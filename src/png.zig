const std = @import("std");

const Error = error {
    WrongMagic,
    WrongDeflateAdler32,
    UnsupportedFilteringType,
    NoPngHeader,
    InconsistentSize,
    UnsupportedInterlaceMethod,
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
    color_type: u8 = 0, // 0, 2, 3, 4, or 6
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
        return @enumFromInt(self.color_type);
    }

    pub fn getPixelSize(self: Self) usize {
        return switch (self.getColorType()) {
            // https://en.wikipedia.org/wiki/PNG#Pixel_format
            ColorType.grayscale => if (self.bit_depth >= 8) (1 * self.bit_depth) / 8
                else @panic("unsupported, packed bit not implemented"),
            ColorType.truecolor => (3 * self.bit_depth) / 8,
            ColorType.indexed => @panic("unsupported, packed bit not implemented"),
            ColorType.grayscaleAlpha => (2 * self.bit_depth) / 8,
            ColorType.truecolorAlpha => (4 * self.bit_depth) / 8,
        };
    }

    pub fn getWidthInBytes(self: Self) usize {
        return self.getWidth() * self.getPixelSize();
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

fn unzip(allocator: std.mem.Allocator, data: std.ArrayList(u8)) ![]u8 {
    var output = std.ArrayList(u8).init(allocator);
    var stream = std.io.fixedBufferStream(data.items[0..]);
    const reader = stream.reader();
    var dcp = std.compress.zlib.decompressor(reader);
    // Decompress
    try dcp.reader().readAllArrayList(&output, 20_000_000_000); // not sure why we need a max value here...
    return try output.toOwnedSlice();
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
fn defiltering(header: IhdrHeader, output: []u8) !usize {
    var line_index: usize = 0; // points at the beginning of the line to defilter
    var write_index: usize = 0; // points at where to write
    const line_width = header.getWidthInBytes();
    while (line_index < output.len) {
        // std.log.debug("line {any}", .{ output[line_index..line_index + line_width + 1] });
        const filtering_type: FilteringType = @enumFromInt(output[line_index]);
        // std.log.debug("defiltering line_index {} filtering_type {}", .{ line_index, filtering_type });
        if (filtering_type == FilteringType.none) {
            std.mem.copyForwards(u8, output[write_index..write_index + line_width], output[line_index + 1..line_index + 1 + line_width]);
            write_index += line_width;
        } else {
            var byte_index = line_index + 1; // points at the next byte to decode
            while ((byte_index - (line_index + 1)) < line_width) {
                // In all the following filtering techniques:
                // A is the equivalent component of the pixel to the left (0 if we are on the first pixel of a line)
                // B is the equivalent component of the pixel to the top (0 if we are on the first line)
                // C is the equivalent component of the pixel to the top left (0 if we are on the first pixel of a line or on the first line)
                const a = if ((write_index % line_width) > 3) output[write_index - 4] else 0;
                const b = if (write_index >= line_width) output[write_index - line_width] else 0;
                const c = if ((write_index >= line_width) and (write_index % line_width) > 3) output[write_index - 4 - line_width] else 0;
                switch (filtering_type) {
                    FilteringType.sub => {
                        output[write_index] = output[byte_index] +% a;
                        byte_index += 1;
                        write_index += 1;
                    },
                    FilteringType.up => {
                        output[write_index] = output[byte_index] +% b;
                        byte_index += 1;
                        write_index += 1;
                    },
                    FilteringType.average => {
                        const a16: u16 = a;
                        const b16: u16 = b;
                        output[write_index] = output[byte_index] +% @as(u8, @intCast((a16 + b16) / 2));
                        byte_index += 1;
                        write_index += 1;
                    },
                    FilteringType.paeth => {
                        const paeth_value: u8 = paeth_predictor(a, b, c);
                        // std.log.debug("byte_index {} write_index {} output[byte_index] {} a {} b {} c {} paeth_value {} -> output[byte_index] +% paeth_value {}", .{ byte_index, write_index, output[byte_index], a, b, c, paeth_value, output[byte_index] +% paeth_value });
                        // std.log.debug("output[byte_index] {} +% paeth_value {} = {}", .{ output[byte_index], paeth_value, output[byte_index] +% paeth_value });
                        output[write_index] = output[byte_index] +% paeth_value;
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

    if (!std.mem.eql(u8, input[0..8], &MAGIC)) {
        return Error.WrongMagic;
    }

    var i = MAGIC.len;
    var end = false;
    var found_header = false;
    while (i < input.len and !end) {
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
                .PLTE => {},
                .IEND => end = true,
                .IDAT => {
                    try idat_data.appendSlice(data);
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

    if (found_header) {
        const uncompressed_data = try unzip(allocator, idat_data);
        const defiltered_size = try defiltering(header, uncompressed_data);

        try idat_data.resize(defiltered_size);
        return PngImageData{
            .header = header,
            .data = uncompressed_data,
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
    var output: [1024]u8 = undefined;
    try load_png(content, &output);
}

test "load_png" {
    std.testing.log_level = .debug;

    var output: [1024]u8 = undefined;
    @memset(&output, 0);
    const no_png_header_error = load_png(&MAGIC, &output);
    try std.testing.expectEqual(Error.NoPngHeader, no_png_header_error);

    @memset(&output, 0);
    _ = try load_png(&MAGIC ++ &[_]u8{
        0x0, 0x0, 0x0, 0xD, // 13 bytes for the header
        'I', 'H', 'D', 'R', // first chunk must be IHDR
        0x0, 0x0, 0x0, 0x0, // width (at 0 for the test to pass)
        0x0, 0x0, 0x0, 0x0, // height (at 0 for the test to pass)
        0x8, // bit depth (1, 2, 4, 8, or 16)
        0x0, // color type (0, 2, 3, 4, or 6)
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
    }, &output);

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

        @memset(&output, 0);
        const header = try load_png(&simple_png, &output);
        const length = header.getWidthInBytes() * header.getHeight();
        try std.testing.expectEqualSlices(u8, &[_]u8{
            // R     G     B     A     R     G     B     A
            0xFF, 0x00, 0x00, 0xFF, 0x00, 0x00, 0xFF, 0xFF,
            0x00, 0xFF, 0x00, 0xFF, 0xFF, 0xFF, 0xFF, 0x00,
        }, output[0..length]);
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

        @memset(&output, 0);
        const header = try load_png(&simple_png, &output);
        const length = header.getWidthInBytes() * header.getHeight();
        try std.testing.expectEqualSlices(u8, &[_]u8{
            // R     G     B     A     R     G     B     A
            0xFF, 0x00, 0x00, 0xFF, 0x00, 0x00, 0xFF, 0xFF,
            0x00, 0xFF, 0x00, 0xFF, 0xFF, 0xFF, 0xFF, 0x00,
        }, output[0..length]);
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

        @memset(&output, 0);
        const header = try load_png(&simple_png, &output);
        const length = header.getWidthInBytes() * header.getHeight();
        try std.testing.expectEqualSlices(u8, &[_]u8{
            // R     G     B     A     R     G     B     A
            0xFF, 0x00, 0x00, 0xFF, 0x00, 0x00, 0xFF, 0xFF,
            0x00, 0xFF, 0x00, 0xFF, 0xFF, 0xFF, 0xFF, 0x00,
        }, output[0..length]);
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

        @memset(&output, 0);
        const header = try load_png(&simple_png, &output);
        const length = header.getWidthInBytes() * header.getHeight();
        try std.testing.expectEqualSlices(u8, &[_]u8{
            // R     G     B     A     R     G     B     A
            0xFF, 0x00, 0x00, 0xFF, 0x00, 0x00, 0xFF, 0xFF,
            0x00, 0xFF, 0x00, 0xFF, 0xFF, 0xFF, 0xFF, 0x00,
        }, output[0..length]);
    }
}

test "adler32" {
    try std.testing.expectEqual(0x11E60398, computeAdler32("Wikipedia"));
}
