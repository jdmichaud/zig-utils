const std = @import("std");

const Error = error {
    WrongMagic,
    WrongDeflateAdler32,
    UnsupportedFilteringType,
    NoPngHeader,
    InconsistentSize,
    UnsupportedInterlaceMethod,
    MissingPalette,
    PaletteFound,
    UnsupportedColorType,
};

//                          P     N     G     \r    \n          \n
const MAGIC = [_]u8{ 0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a };

const ChunkHeader = extern struct {
    length_big_endian: u32, // only include data, must add 12 for the whole chunk
    chunk_type: [4]u8,

    const Self = @This();
    fn getLength(self: Self) u32 {
        return @byteSwap(self.length_big_endian);
    }

    fn getType(self: Self) ?ChunkType {
        return std.meta.stringToEnum(ChunkType, &self.chunk_type);
    }
};

const ChunkType = enum {
    IHDR,
    PLTE,
    IDAT,
    IEND,
    // OPTIONAL
    bKGD,
    cHRM,
    cICP,
    dSIG,
    eXIf,
    gAMA,
    hIST,
    iCCP,
    iTXt,
    pHYs,
    sBIT,
    sPLT,
    sRGB,
    sTER,
    tEXt,
    tIME,
    tRNS,
    zTXt,
};

const IhdrHeader = packed struct {
    width: u32 = 0,
    height: u32 = 0,
    bit_depth: u8 = 0, // 1, 2, 4, 8, or 16
    color_type_tag: u8 = 0, // 0, 2, 3, 4, or 6
    compression_method: u8 = 0,
    filter_method: u8 = 0,
    interlace_method: u8 = 0, // 0 "no interlace" or 1 "Adam7 interlace"

    const Self = @This();

    pub fn getWidth(self: Self) u32 {
        return @byteSwap(self.width);
    }

    pub fn getHeight(self: Self) u32 {
        return @byteSwap(self.height);
    }

    pub fn getColorType(self: Self) ColorType {
        return @enumFromInt(self.color_type_tag);
    }

    pub fn getWidthInBytes(self: Self) usize {
        return switch (self.getColorType()) {
            // https://en.wikipedia.org/wiki/PNG#Pixel_format
            ColorType.grayscale => if (self.bit_depth >= 8) self.getWidth() * (1 * self.bit_depth) / 8
                else @panic("unsupported, packed bit not implemented"),
            ColorType.truecolor => self.getWidth() * (3 * self.bit_depth) / 8,
            ColorType.indexed => ((self.getWidth() * self.bit_depth) + 7) / 8,
            ColorType.grayscaleAlpha => self.getWidth() * (2 * self.bit_depth) / 8,
            ColorType.truecolorAlpha => self.getWidth() * (4 * self.bit_depth) / 8,
        };
    }
};

pub const ColorType = enum(u8) {
    grayscale = 0,
    truecolor = 2,
    indexed = 3,
    grayscaleAlpha = 4,
    truecolorAlpha = 6,
};

// https://en.wikipedia.org/wiki/Adler-32#Calculation
fn computeAdler32(data: []const u8) u32 {
    const modulo = 65521;
    var A: u32 = 1;
    var B: u32 = 0;

    for (data) |byte| {
        A = (A + byte) % modulo;
        B = (A + B) % modulo;
    }

    return B << 16 | A;
}

const HuffmanLengthRange = struct {
    symbolStart: u16,   // first literal/length symbol in the range
    symbolEnd: u16,     // last literal/length symbol in the range
    codeLength: u8,     // length of the Huffman code in bits
};

// Fixed Huffman codes as described in RFC 1951
// Lit Value    Bits        Codes
// ---------    ----        -----
//  0 - 143     8           00110000 through
//                          10111111
// 144 - 255     9          110010000 through
//                          111111111
// 256 - 279     7          0000000 through
//                          0010111
// 280 - 287     8          11000000 through
//                          11000111
const FixedHuffmanTable: []HuffmanLengthRange = &[_]HuffmanLengthRange{
    .{ .symbolStart = 0,   .symbolEnd = 143, .codeLength = 8 },
    .{ .symbolStart = 144, .symbolEnd = 255, .codeLength = 9 },
    .{ .symbolStart = 256, .symbolEnd = 279, .codeLength = 7 },
    .{ .symbolStart = 280, .symbolEnd = 287, .codeLength = 8 },
};

fn unzip(allocator: std.mem.Allocator, data: std.ArrayList(u8)) !std.ArrayList(u8) {
    var output = std.ArrayList(u8).init(allocator);
    var stream = std.io.fixedBufferStream(data.items[0..]);
    const reader = stream.reader();
    var dcp = std.compress.zlib.decompressor(reader);
    // Decompress
    try dcp.reader().readAllArrayList(&output, 20_000_000_000); // not sure why we need a max value here...
    return output;
}

const FilteringType = enum(u8) {
    none = 0,
    sub = 1,
    up = 2,
    average = 3,
    paeth = 4,
};

fn paeth_predictor(a: u8, b: u8, c: u8) u8 {
    const v_8: @Vector(3, u8) = .{ a, b, c };
    const v_16: @Vector(3, i16) = @intCast(v_8);
    const p: @Vector(3, i16) = @splat(v_16[0] + v_16[1] - v_16[2]); // initial estimate
    const distance: @Vector(3, u16) = @abs(p - v_16); // distances to a, b, c

    if (distance[0] <= distance[1] and distance[0] <= distance[2]) {
        return a;
    } else if (distance[1] <= distance[2]) {
        return b;
    } else {
        return c;
    }
}

// PNG defiltering as described here: https://en.wikipedia.org/wiki/PNG#Filtering
// Remember that filtering is applied on bytes not pixels.
fn defiltering(header: IhdrHeader, output: std.ArrayList(u8)) !usize {
    var line_index: usize = 0; // points at the beginning of the line to defilter
    var write_index: usize = 0; // points at where to write
    const line_width = header.getWidthInBytes();
    while (line_index < output.items.len) {
        const filtering_type: FilteringType = @enumFromInt(output.items[line_index]);
        if (filtering_type == FilteringType.none) {
            std.mem.copyForwards(u8, output.items[write_index..write_index + line_width], output.items[line_index + 1..line_index + 1 + line_width]);
            write_index += line_width;
        } else {
            var byte_index = line_index + 1; // points at the next byte to decode
            while ((byte_index - (line_index + 1)) < line_width) {
                // In all the following filtering techniques:
                // A is the equivalent component of the pixel to the left (0 if we are on the first pixel of a line)
                // B is the equivalent component of the pixel to the top (0 if we are on the first line)
                // C is the equivalent component of the pixel to the top left (0 if we are on the first pixel of a line or on the first line)
                const a = if ((write_index % line_width) > 3) output.items[write_index - 4] else 0;
                const b = if (write_index >= line_width) output.items[write_index - line_width] else 0;
                const c = if ((write_index >= line_width) and (write_index % line_width) > 3) output.items[write_index - 4 - line_width] else 0;
                switch (filtering_type) {
                    FilteringType.sub => {
                        output.items[write_index] = output.items[byte_index] +% a;
                        byte_index += 1;
                        write_index += 1;
                    },
                    FilteringType.up => {
                        output.items[write_index] = output.items[byte_index] +% b;
                        byte_index += 1;
                        write_index += 1;
                    },
                    FilteringType.average => {
                        const a16: u16 = a;
                        const b16: u16 = b;
                        output.items[write_index] = output.items[byte_index] +% @as(u8, @intCast((a16 + b16) / 2));
                        byte_index += 1;
                        write_index += 1;
                    },
                    FilteringType.paeth => {
                        const paeth_value: u8 = paeth_predictor(a, b, c);
                        output.items[write_index] = output.items[byte_index] +% paeth_value;
                        byte_index += 1;
                        write_index += 1;
                    },
                    else => {
                        std.log.warn("unsupported filter type {}", .{ filtering_type });
                        return Error.UnsupportedFilteringType;
                    },
                }
            }
        }
        line_index += line_width + 1;
    }
    // std.log.debug("line_width {} header.getHeight() {}", .{ line_width, header.getHeight() });
    // std.log.debug("write_index {} line_width * header.getHeight() {}", .{ write_index, line_width * header.getHeight() });
    if (write_index != line_width * header.getHeight()) {
        return Error.InconsistentSize;
    }

    return write_index;
}

fn indexed_to_rgba(allocator: std.mem.Allocator, header: IhdrHeader, palette: []const u8,
    transparency: []const u8, uncompressed_data: []const u8) !std.ArrayList(u8) {
    var converted = std.ArrayList(u8).init(allocator);

    // We want to extract u2 values
    const num_values = header.getWidth() * header.getHeight();

    var bit_offset: usize = 0;
    var i: usize = 0;
    // The amount of bit to pad at the end of each rows
    const line_padding = (header.getWidth() * header.bit_depth) % 8;
    const mask: u8 = switch (header.bit_depth) {
        1 => 0b10000000,
        2 => 0b11000000,
        4 => 0b11110000,
        8 => 0b11111111,
        else => unreachable,
    };
    while (i < num_values) : (i += 1) {
        const index = bit_offset / 8;
        const mask_shift: u3 = @intCast(bit_offset % 8);
        const result_shift: u3 = @intCast(8 - header.bit_depth - mask_shift);
        const val: usize = (uncompressed_data[index] & (mask >> mask_shift)) >> result_shift;
        // std.log.debug("bit_offset {} line_padding {} index {} mask >> mask_shift {b} result_shift {} => val {}", .{
        //     bit_offset, line_padding, index, mask >> mask_shift, result_shift, val
        // });

        try converted.append(palette[val * 3]);
        try converted.append(palette[val * 3 + 1]);
        try converted.append(palette[val * 3 + 2]);
        try converted.append(transparency[val]);

        if ((i + 1) % header.getWidth() == 0) {
            // The beginning of a row ia always byte aligned so the last byte of
            // a row can be padded
            bit_offset += header.bit_depth + line_padding;
        } else {
            bit_offset += header.bit_depth;
        }
    }

    return converted;
}

pub const PngImageData = struct {
    header: IhdrHeader,
    data: []u8,

    const Self = @This();

    pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
        allocator.free(self.data);
    }
};

// PNG → IDAT → zlib → DEFLATE
pub fn load_png(allocator: std.mem.Allocator, input: []const u8) !PngImageData {
    var idat_data = std.ArrayList(u8).init(allocator);
    defer idat_data.deinit();
    var header: IhdrHeader = .{};
    var palette: [256 * 3]u8 = undefined;
    var transparency: [256]u8 = [_]u8{ 0xFF } ** 256;

    if (!std.mem.eql(u8, input[0..8], &MAGIC)) {
        return Error.WrongMagic;
    }

    var i = MAGIC.len;
    var end = false;
    var found_header = false;
    var found_palette = false;
    while (i < input.len and !end) {
        // FIXME check CRC
        const chunk_header: *align(1) const ChunkHeader = @alignCast(@ptrCast(input[i..i + @sizeOf(ChunkHeader)].ptr));
        std.log.debug("chunk_type {s} chunk_length 0x{x}", .{ chunk_header.chunk_type, chunk_header.getLength() });
        const data = input[i + 8..i + 8 + chunk_header.getLength()];
        if (chunk_header.getType()) |chunk_type| {
            switch (chunk_type) {
                .IHDR => {
                    found_header = true;
                    const hdr_bytes: [*]u8 = @ptrCast(&header);
                    @memcpy(hdr_bytes, data[0..13]);
                    std.log.debug("IHDR width {} height {} bit depth {} color type {} compression method {}", .{
                        header.getWidth(),
                        header.getHeight(),
                        header.bit_depth,
                        header.getColorType(),
                        header.compression_method,
                    });
                    if (header.interlace_method != 0) {
                        return Error.UnsupportedInterlaceMethod;
                    }
                },
                .PLTE => {
                    found_palette = true;
                    @memcpy(palette[0..data.len], data);
                },
                .IEND => end = true,
                .IDAT => {
                    try idat_data.appendSlice(data);
                },
                .tRNS => {
                    @memcpy(transparency[0..data.len], data);
                },
                else => {
                    std.log.debug("chunk {any} ignored", .{ chunk_header.getType() });
                }
            }
        } else {
            std.log.debug("unknown chunk type {s}", .{ chunk_header.chunk_type });
        }
        i += chunk_header.getLength() + 4 + 4 + 4; // data + byte length + type + CRC
    }

    if (!end) {
        std.log.warn("png stream finished without a END chunk", .{});
    }

    if (!found_palette and header.getColorType() == ColorType.indexed) {
        return Error.MissingPalette;
    }

    if (found_header) {
        var uncompressed_data = try unzip(allocator, idat_data);
        const defiltered_size = try defiltering(header, uncompressed_data);
        try uncompressed_data.resize(defiltered_size);

        const data = brk: switch (header.getColorType()) {
            ColorType.truecolorAlpha => break :brk try uncompressed_data.toOwnedSlice(),
            ColorType.indexed => {
                var converted_data = try indexed_to_rgba(allocator, header, &palette, &transparency, uncompressed_data.items);
                uncompressed_data.deinit();
                break :brk try converted_data.toOwnedSlice();
            },
            else => {
                std.log.err("unsupported color type {}", .{ header.getColorType() });
                return Error.UnsupportedColorType;
            },
        };

        return PngImageData{
            .header = header,
            .data = data,
        };
    }

    return Error.NoPngHeader;
}

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const args = try std.process.argsAlloc(allocator);
    defer std.process.argsFree(allocator, args);

    const content = try std.fs.cwd().readFileAlloc(allocator, args[1], std.math.maxInt(usize));
    defer allocator.free(content);
    const png = try load_png(allocator, content);
    png.deinit(allocator);
}

test "load_png" {
    std.testing.log_level = .debug;

    const no_png_header_error = load_png(std.testing.allocator, &MAGIC);
    try std.testing.expectEqual(Error.NoPngHeader, no_png_header_error);

    {
        var png = try load_png(std.testing.allocator, &MAGIC ++ &[_]u8{
            0x0, 0x0, 0x0, 0xD, // 13 bytes for the header
            'I', 'H', 'D', 'R', // first chunk must be IHDR
            0x0, 0x0, 0x0, 0x0, // width (at 0 for the test to pass)
            0x0, 0x0, 0x0, 0x0, // height (at 0 for the test to pass)
            0x8, // bit depth (1, 2, 4, 8, or 16)
            0x6, // color type (0, 2, 3, 4, or 6)
            0x0, // compression method (always 0)
            0x0, // filter method (always 0)
            0x0, // interlace method (1 byte, values 0 "no interlace" or 1 "Adam7 interlace")
            0x0, 0x0, 0x0, 0x0, // CRC

            0x0, 0x0, 0x0, 0x8, // size
            'I', 'D', 'A', 'T',
            // Simplest zlib payload possible
            // 78 → CMF byte
            // - CM = 8 → DEFLATE
            // - CINFO = 7 → 32 KB window
            // 9C → FLG byte
            // - Compression level = default
            // - Check bits make (CMF*256 + FLG) % 31 == 0
            // 0x03 means:
            // - BFINAL = 1 → last block
            // - BTYPE = 01 → fixed Huffman
            // 00 → end-of-block symbol (encoded by fixed Huffman)
            // 00 00 00 01 → Adlter32 of an empty block is 1
            0x78, 0x9C, 0x03, 0x00, 0x00, 0x00, 0x00, 0x01,
            0x0, 0x0, 0x0, 0x0, // CRC

            0x0, 0x0, 0x0, 0x0, // size of the data (does include chunk type and CRC)
            'I', 'E', 'N', 'D', // last chunk must be IEND
            0x0, 0x0, 0x0, 0x0, // CRC
        });
        png.deinit(std.testing.allocator);
    }

    {
        // A simple 2x2 png layout as follow:
        // red, blue
        // green, transparent white
        // Using a paeth filter for the second row
        const simple_png = [_]u8{
          0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A, // magic
          0x00, 0x00, 0x00, 0x0D, // 13 bytes for the header
          0x49, 0x48, 0x44, 0x52, // first chunk must be IHDR
          0x00, 0x00, 0x00, 0x02, // width
          0x00, 0x00, 0x00, 0x02, // height
          0x08, // bit depth (1, 2, 4, 8, or 16)
          0x06, // color type (0, 2, 3, 4, or 6)
          0x00, // compression method (always 0)
          0x00, // filter method (always 0)
          0x00, // interlace method (1 byte, values 0 "no interlace" or 1 "Adam7 interlace")
          0x72, 0xB6, 0x0D, 0x24, // CRC
          0x00, 0x00, 0x00, 0x16, // size: 22
          0x49, 0x44, 0x41, 0x54, // "IDAT"
          0x78, 0x9C, 0x63, 0xF8, 0xCF, 0xC0, 0xF0, 0x9F, 0x81, 0xE1, 0xFF, 0x7F, 0x16, 0x46, 0x08, 0x8B, 0x11, 0x00, 0x3F, 0x00, 0x06, 0x01,
          0x3C, 0x23, 0xB5, 0xE9, // CRC
          0x00, 0x00, 0x00, 0x00, // size: 0
          0x49, 0x45, 0x4E, 0x44, // "IEND"
          0xAE, 0x42, 0x60, 0x82, // CRC
        };

        var png = try load_png(std.testing.allocator, &simple_png);
        try std.testing.expectEqualSlices(u8, &[_]u8{
            // R     G     B     A     R     G     B     A
            0xFF, 0x00, 0x00, 0xFF, 0x00, 0x00, 0xFF, 0xFF,
            0x00, 0xFF, 0x00, 0xFF, 0xFF, 0xFF, 0xFF, 0x00,
        }, png.data);
        png.deinit(std.testing.allocator);
    }

    {
        // A simple 2x2 png layout as follow:
        // red, blue
        // green, transparent white
        // Using a sub filter for the second row
        const simple_png = [_]u8{
          0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A, // magic
          0x00, 0x00, 0x00, 0x0D, // 13 bytes for the header
          0x49, 0x48, 0x44, 0x52, // first chunk must be IHDR
          0x00, 0x00, 0x00, 0x02, // width
          0x00, 0x00, 0x00, 0x02, // height
          0x08, // bit depth (1, 2, 4, 8, or 16)
          0x06, // color type (0, 2, 3, 4, or 6)
          0x00, // compression method (always 0)
          0x00, // filter method (always 0)
          0x00, // interlace method (1 byte, values 0 "no interlace" or 1 "Adam7 interlace")
          0x72, 0xB6, 0x0D, 0x24, // CRC
          0x00, 0x00, 0x00, 0x15, // size: 21
          0x49, 0x44, 0x41, 0x54, // "IDAT"
          0x78, 0x9C, 0x63, 0xF8, 0xCF, 0xC0, 0x00, 0x42, 0xFF, 0x19, 0x81, 0xD4, 0x7F, 0x20, 0x62, 0x04, 0x00, 0x45, 0xD6, 0x07, 0xFB,
          0x76, 0x4D, 0x98, 0x2C, // CRC
          0x00, 0x00, 0x00, 0x00, // size: 0
          0x49, 0x45, 0x4E, 0x44, // "IEND"
          0xAE, 0x42, 0x60, 0x82, // CRC
        };

        var png = try load_png(std.testing.allocator, &simple_png);
        try std.testing.expectEqualSlices(u8, &[_]u8{
            // R     G     B     A     R     G     B     A
            0xFF, 0x00, 0x00, 0xFF, 0x00, 0x00, 0xFF, 0xFF,
            0x00, 0xFF, 0x00, 0xFF, 0xFF, 0xFF, 0xFF, 0x00,
        }, png.data);
        png.deinit(std.testing.allocator);
    }

    {
        // A simple 2x2 png layout as follow:
        // red, blue
        // green, transparent white
        // Using an average filter for the second row
        const simple_png = [_]u8{
          0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A, // magic
          0x00, 0x00, 0x00, 0x0D, // 13 bytes for the header
          0x49, 0x48, 0x44, 0x52, // first chunk must be IHDR
          0x00, 0x00, 0x00, 0x02, // width
          0x00, 0x00, 0x00, 0x02, // height
          0x08, // bit depth (1, 2, 4, 8, or 16)
          0x06, // color type (0, 2, 3, 4, or 6)
          0x00, // compression method (always 0)
          0x00, // filter method (always 0)
          0x00, // interlace method (1 byte, values 0 "no interlace" or 1 "Adam7 interlace")
          0x72, 0xB6, 0x0D, 0x24, // CRC
          0x00, 0x00, 0x00, 0x18, // Size: 24 bytes
          0x49, 0x44, 0x41, 0x54, // "IDAT"
          0x78, 0x9C, 0x63, 0xF8, 0xCF, 0xC0, 0x00, 0x42, 0xFF, 0x99, 0x1B, 0xFF, 0x33, 0x34, 0xFC,
          0x6F, 0x68, 0x60, 0x04, 0x00, 0x47, 0xF7, 0x08, 0x00,
          0x24, 0x9C, 0x6E, 0x42, // CRC
          0x00, 0x00, 0x00, 0x00, // size: 0
          0x49, 0x45, 0x4E, 0x44, // "IEND"
          0xAE, 0x42, 0x60, 0x82, // CRC
        };

        var png = try load_png(std.testing.allocator, &simple_png);
        try std.testing.expectEqualSlices(u8, &[_]u8{
            // R     G     B     A     R     G     B     A
            0xFF, 0x00, 0x00, 0xFF, 0x00, 0x00, 0xFF, 0xFF,
            0x00, 0xFF, 0x00, 0xFF, 0xFF, 0xFF, 0xFF, 0x00,
        }, png.data);
        png.deinit(std.testing.allocator);
    }

    {
        // A simple 2x2 png layout as follow:
        // red, blue
        // green, transparent white
        // Using a "up" filter for the second row
        const simple_png = [_]u8{
          0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A, // magic
          0x00, 0x00, 0x00, 0x0D, // 13 bytes for the header
          0x49, 0x48, 0x44, 0x52, // first chunk must be IHDR
          0x00, 0x00, 0x00, 0x02, // width
          0x00, 0x00, 0x00, 0x02, // height
          0x08, // bit depth (1, 2, 4, 8, or 16)
          0x06, // color type (0, 2, 3, 4, or 6)
          0x00, // compression method (always 0)
          0x00, // filter method (always 0)
          0x00, // interlace method (1 byte, values 0 "no interlace" or 1 "Adam7 interlace")
          0x72, 0xB6, 0x0D, 0x24, // CRC
          0x00, 0x00, 0x00, 0x14, // Size: 20 bytes
          0x49, 0x44, 0x41, 0x54, // "IDAT"
          0x78, 0x9C, 0x63, 0xF8, 0xCF, 0xC0, 0x00, 0x42,  0xFF, 0x99, 0x18, 0xC1, 0x14, 0x03, 0x23,
          0x00, 0x41, 0xEB, 0x06, 0xFE,
          0xC5, 0xFC, 0x7D, 0x6F, // CRC
          0x00, 0x00, 0x00, 0x00, // size: 0
          0x49, 0x45, 0x4E, 0x44, // "IEND"
          0xAE, 0x42, 0x60, 0x82, // CRC
        };

        var png = try load_png(std.testing.allocator, &simple_png);
        try std.testing.expectEqualSlices(u8, &[_]u8{
            // R     G     B     A     R     G     B     A
            0xFF, 0x00, 0x00, 0xFF, 0x00, 0x00, 0xFF, 0xFF,
            0x00, 0xFF, 0x00, 0xFF, 0xFF, 0xFF, 0xFF, 0x00,
        }, png.data);
        png.deinit(std.testing.allocator);
    }

    {
        // A simple 2x2 png layout as follow:
        // red, blue
        // green, transparent white
        // Using a indexed pixel format with a color palette 2 bits of depth
        const simple_png = [_]u8{
          0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A, // magic
          0x00, 0x00, 0x00, 0x0D, // 13 bytes for the header
          0x49, 0x48, 0x44, 0x52, // first chunk must be IHDR
          0x00, 0x00, 0x00, 0x02, // width
          0x00, 0x00, 0x00, 0x02, // height
          0x02, // bit depth (1, 2, 4, 8, or 16)
          0x03, // indexed
          0x00, // compression method (always 0)
          0x00, // filter method (always 0)
          0x00, // interlace method (1 byte, values 0 "no interlace" or 1 "Adam7 interlace")
          0x72, 0xB6, 0x0D, 0x24, // CRC
          // =====================
          // PLTE (4 colors)
          // =====================
          0x00, 0x00, 0x00, 0x0C,             // length = 12
          0x50, 0x4C, 0x54, 0x45,             // "PLTE"
              0xFF, 0x00, 0x00,               // red
              0x00, 0x00, 0xFF,               // blue
              0x00, 0xFF, 0x00,               // green
              0xFF, 0xFF, 0xFF,               // white (transparent via tRNS)
          0xF0, 0x50, 0x0B, 0x0F,             // CRC
          // =====================
          // tRNS (alpha for 4 entries)
          // =====================
          0x00, 0x00, 0x00, 0x04,             // length = 4
          0x74, 0x52, 0x4E, 0x53,             // "tRNS"
              0xFF, 0xFF, 0xFF, 0x00,         // last entry fully transparent
          0x6B, 0xFA, 0xB0, 0x26,             // CRC

          // =====================
          // IDAT (Image Data)
          // =====================
          0x00, 0x00, 0x00, 0x0F, // Length: 15 bytes
          0x49, 0x44, 0x41, 0x54, // "IDAT"

          // ZLIB STREAM
          0x78, 0x01,             // Zlib Header
          0x01,                   // Block Header (Final=1, Type=00 Stored)
          0x04, 0x00,             // Len = 4 bytes
          0xFB, 0xFF,             // NLen (One's complement)

          // RAW DATA
          0x00,                   // Filter Row 1 (None)
          0x10,                   // Row 1: Red(00) Blue(01) -> 0x10
          0x00,                   // Filter Row 2 (None)
          0xB0,                   // Row 2: Green(10) Trans(11) -> 0xB0

          // ADLER-32 CHECKSUM
          // Calculated on "00 10 00 B0"
          0x00, 0xE4, 0x00, 0xC1,

          // IDAT CRC-32
          // Calculated on "IDAT" + Zlib Stream
          0x5C, 0x36, 0x0B, 0x92,

          // =====================
          // IEND
          // =====================
          0x00, 0x00, 0x00, 0x00,
          0x49, 0x45, 0x4E, 0x44,
          0xAE, 0x42, 0x60, 0x82,
        };
        var png = try load_png(std.testing.allocator, &simple_png);
        try std.testing.expectEqualSlices(u8, &[_]u8{
            // R     G     B     A     R     G     B     A
            0xFF, 0x00, 0x00, 0xFF, 0x00, 0x00, 0xFF, 0xFF,
            0x00, 0xFF, 0x00, 0xFF, 0xFF, 0xFF, 0xFF, 0x00,
        }, png.data);
        png.deinit(std.testing.allocator);
    }

    {
        // A simple 2x2 png layout as follow:
        // red, blue
        // green, transparent white
        // Using a indexed pixel format with a color palette 8 bits of depth
        const simple_png = [_]u8{
            0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A, // magic
            0x00, 0x00, 0x00, 0x0D, // 13 bytes for the header
            0x49, 0x48, 0x44, 0x52, // IHDR
            0x00, 0x00, 0x00, 0x02, // width
            0x00, 0x00, 0x00, 0x02, // height
            0x08, // bit depth 8
            0x03, // indexed
            0x00, // compression
            0x00, // filter
            0x00, // interlace
            0x29, 0xF7, 0x3E, 0x21,

            // =====================
            // PLTE (4 colors)
            // =====================
            0x00, 0x00, 0x00, 0x0C,
            0x50, 0x4C, 0x54, 0x45,
                0xFF, 0x00, 0x00, // 0: Red
                0x00, 0x00, 0xFF, // 1: Blue
                0x00, 0xFF, 0x00, // 2: Green
                0xFF, 0xFF, 0xFF, // 3: White/Trans
            0xF0, 0x50, 0x0B, 0x0F,

            // =====================
            // tRNS (alpha)
            // =====================
            0x00, 0x00, 0x00, 0x04,
            0x74, 0x52, 0x4E, 0x53,
            0xFF, 0xFF, 0xFF, 0x00,
            0x6B, 0xFA, 0xB0, 0x26,

            // =====================
            // IDAT (Image Data)
            // =====================
            // Previous Length: 15. New Length: 17
            // (Zlib header 2 + Block header 5 + Raw data 6 + Adler 4)
            0x00, 0x00, 0x00, 0x11, // <--- CHANGED: Length
            0x49, 0x44, 0x41, 0x54, // "IDAT"

            // ZLIB STREAM
            0x78, 0x01,             // Zlib Header
            0x01,                   // Block Header (Final=1, Type=00 Stored)
            0x06, 0x00,             // <--- CHANGED: Len = 6 bytes (2 rows * 3 bytes)
            0xF9, 0xFF,             // <--- CHANGED: NLen (One's complement of 0x0006)

            // RAW DATA
            0x00,                   // Filter Row 1 (None)
            0x00, 0x01,             // Row 1: Index 00, Index 01
            0x00,                   // Filter Row 2 (None)
            0x02, 0x03,             // Row 2: Index 02, Index 03

            // ADLER-32 CHECKSUM
            0x00, 0x11, 0x00, 0x07,

            // IDAT CRC-32
            0xC0, 0x99, 0x76, 0x69,

            // =====================
            // IEND - Unchanged
            // =====================
            0x00, 0x00, 0x00, 0x00,
            0x49, 0x45, 0x4E, 0x44,
            0xAE, 0x42, 0x60, 0x82,
        };
        var png = try load_png(std.testing.allocator, &simple_png);
        try std.testing.expectEqualSlices(u8, &[_]u8{
            // R     G     B     A     R     G     B     A
            0xFF, 0x00, 0x00, 0xFF, 0x00, 0x00, 0xFF, 0xFF,
            0x00, 0xFF, 0x00, 0xFF, 0xFF, 0xFF, 0xFF, 0x00,
        }, png.data);
        png.deinit(std.testing.allocator);
    }
}

test "adler32" {
    try std.testing.expectEqual(0x11E60398, computeAdler32("Wikipedia"));
}
