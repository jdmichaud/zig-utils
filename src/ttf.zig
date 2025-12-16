// https://developer.apple.com/fonts/TrueType-Reference-Manual/RM06/Chap6.html
const std = @import("std");

/// Reads a T in big endian from b at index and advance the index
fn readT(comptime T: type, b: []const u8, index: *usize) T {
  const size = @sizeOf(T);
  const v = std.mem.readInt(T, b[index.*..][0..size], .big);
  index.* += size;
  return v;
}

const OffsetTable = struct {
  scaler_type: u32,     // A tag to indicate the OFA scaler to be used to rasterize this font;
                        // see the note on the scaler type below for more information.
  num_tables: u16,      // number of tables
  search_range: u16,    // (maximum power of 2 <= numTables)*16
  entry_selector: u16,  // log2(maximum power of 2 <= numTables)
  range_shift: u16,     // numTables*16-searchRange

  fn fromBytes(bytes: []const u8, off: *usize) OffsetTable {
    return .{
      .scaler_type = readT(u32, bytes, off),
      .num_tables = readT(u16, bytes, off),
      .search_range = readT(u16, bytes, off),
      .entry_selector = readT(u16, bytes, off),
      .range_shift = readT(u16, bytes, off),
    };
  }

};

const TableRecord = struct {
  tag: [4]u8,
  checksum: u32,
  offset: u32,
  length: u32,

  fn fromBytes(bytes: []const u8, off: *usize) TableRecord {
    var tag: [4]u8 = undefined;
    @memcpy(&tag, bytes[off.*..][0..4]);
    off.* += 4;

    return .{
      .tag = tag,
      .checksum = readT(u32, bytes, off),
      .offset = readT(u32, bytes, off),
      .length = readT(u32, bytes, off),
    };
  }
};

/// Helper type for 16.16 fixed-point values
pub const Fixed = struct {
  raw: i32,

  pub fn fromParts(integer: i16, fraction: u16) Fixed {
    return .{ .raw = (@as(i32, integer) << 16) | fraction };
  }

  pub fn toFloat(self: Fixed) f32 {
    return @as(f32, self.raw) / 65536.0;
  }
};

pub const HeadTable = struct {
  /// 16.16 fixed-point number
  /// Always 0x00010000 if (version 1.0)
  version: Fixed,
  /// 16.16 fixed-point number
  /// set by font manufacturer
  font_revision: Fixed,

  check_sum_adjustment: u32,
  /// set to 0x5F0F3CF5
  magic_number: u32,
  /// bit 0 - y value of 0 specifies baseline
  /// bit 1 - x position of left most black bit is LSB
  /// bit 2 - scaled point size and actual point size will differ (i.e. 24 point glyph differs from 12 point glyph scaled by factor of 2)
  /// bit 3 - use integer scaling instead of fractional
  /// bit 4 - (used by the Microsoft implementation of the TrueType scaler)
  /// bit 5 - This bit should be set in fonts that are intended to e laid out vertically, and in which the glyphs have been drawn such that an x-coordinate of 0 corresponds to the desired vertical baseline.
  /// bit 6 - This bit must be set to zero.
  /// bit 7 - This bit should be set if the font requires layout for correct linguistic rendering (e.g. Arabic fonts).
  /// bit 8 - This bit should be set for an AAT font which has one or more metamorphosis effects designated as happening by default.
  /// bit 9 - This bit should be set if the font contains any strong right-to-left glyphs.
  /// bit 10 - This bit should be set if the font contains Indic-style rearrangement effects.
  /// bits 11-13 - Defined by Adobe.
  /// bit 14 - This bit should be set if the glyphs in the font are simply generic symbols for code point ranges, such as for a last resort font.
  flags: u16,
  /// range from 64 to 16384
  units_per_em: u16,
  /// Seconds since 1904-01-01 (Mac epoch)
  created: i64,
  /// Seconds since 1904-01-01 (Mac epoch)
  modified: i64,
  x_min: i16, // for all glyph bounding boxes
  y_min: i16, // for all glyph bounding boxes
  x_max: i16, // for all glyph bounding boxes
  y_max: i16, // for all glyph bounding boxes
  /// bit 0 bold
  /// bit 1 italic
  /// bit 2 underline
  /// bit 3 outline
  /// bit 4 shadow
  /// bit 5 condensed (narrow)
  /// bit 6 extended
  mac_style: u16,
  /// smallest readable size in pixels
  lowest_rec_ppem: u16,
  ///  0 Mixed directional glyphs
  ///  1 Only strongly left to right glyphs
  ///  2 Like 1 but also contains neutrals
  /// -1 Only strongly right to left glyphs
  /// -2 Like -1 but also contains neutrals
  font_direction_hint: i16,
  /// 0 = short offsets, 1 = long offsets
  index_to_loc_format: i16,
  /// Must be 0
  glyph_data_format: i16,

  // Size of the head table
  pub const byte_len: usize = 54;
  /// Parse a `head` table from a big-endian byte slice
  pub fn init(bytes: []const u8) !HeadTable {
    if (bytes.len < byte_len) {
      return error.UnexpectedEndOfTable;
    }

    var off: usize = 0;

    return HeadTable{
      .version = .{ .raw = readT(i32, bytes, &off) },
      .font_revision = .{ .raw = readT(i32, bytes, &off) },

      .check_sum_adjustment = readT(u32, bytes, &off),
      .magic_number = readT(u32, bytes, &off),

      .flags = readT(u16, bytes, &off),
      .units_per_em = readT(u16, bytes, &off),

      .created = readT(i64, bytes, &off),
      .modified = readT(i64, bytes, &off),

      .x_min = readT(i16, bytes, &off),
      .y_min = readT(i16, bytes, &off),
      .x_max = readT(i16, bytes, &off),
      .y_max = readT(i16, bytes, &off),

      .mac_style = readT(u16, bytes, &off),
      .lowest_rec_ppem = readT(u16, bytes, &off),

      .font_direction_hint = readT(i16, bytes, &off),
      .index_to_loc_format = readT(i16, bytes, &off),
      .glyph_data_format = readT(i16, bytes, &off),
    };
  }
};

const CmapHeader = struct {
  /// Version number (Set to zero)
  version: u16,
  /// Number of encoding subtables
  numberSubtables: u16,
};

const PlatformID = enum(u16) {
  /// Indicates Unicode version
  unicode = 0,
  /// Script Manager code
  macintosh = 1,
  /// Do not use
  reserved = 2,
  /// Microsoft encoding
  microsoft = 3,
};

const CmapRecord = struct {
  /// Platform identifier
  platformID: PlatformID,
  /// Platform-specific encoding identifier
  platformSpecificID: u16,
  /// Offset of the mapping table
  offset: u32,
};

const CmapSubtable = struct {
  format: u16,
  /// length in bytes including this structure
  length: u16,
  /// Language code
  language: u16,
  /// An array that maps character codes to glyph index values
  glyph_index: []const u8,
};

pub const CmapTable = struct {
  header: CmapHeader,
  records: []CmapRecord,
  tables: []CmapSubtable,

  const Self = @This();
  pub fn init(allocator: std.mem.Allocator, bytes: []const u8) !CmapTable {
    var index: usize = 0;

    // 1. Read the Header
    const header = CmapHeader{
      .version = readT(u16, bytes, &index),
      .numberSubtables = readT(u16, bytes, &index),
    };

    // 2. Allocate and Read Encoding Records
    const records = try allocator.alloc(CmapRecord, header.numberSubtables);
    errdefer allocator.free(records);
    for (records) |*rec| {
      rec.platformID = @enumFromInt(readT(u16, bytes, &index));
      rec.platformSpecificID = readT(u16, bytes, &index);
      rec.offset = readT(u32, bytes, &index);
    }

    // 3. Read Subtables
    // Note: The records contain offsets relative to the start of the Cmap Table.
    // We must jump to those offsets to read the actual subtable data.
    const tables = try allocator.alloc(CmapSubtable, header.numberSubtables);
    errdefer allocator.free(tables);

    for (records, 0..) |rec, i| {
      // Create a cursor at the offset for this specific subtable
      var sub_cursor = @as(usize, rec.offset);

      // Check bounds
      if (sub_cursor >= bytes.len) return error.InvalidOffset;

      // Peek at the format (first u16) to determine how to read the length
      // We can't use readT here initially because we don't want to advance sub_cursor yet
      const format = readT(u16, bytes, &sub_cursor);

      var length: u32 = 0;
      var language: u16 = 0;
      var header_size: usize = 0;

      // Handle standard formats vs Format 12 (High capacity)
      if (format == 12 or format == 13) {
        // Format 12 Header: format(u16) + reserved(u16) + length(u32) + language(u32) + ...
        _ = readT(u16, bytes, &sub_cursor); // skip reserved
        length = readT(u32, bytes, &sub_cursor);
        language = @truncate(readT(u32, bytes, &sub_cursor)); // Truncate u32 lang to match your struct
        header_size = 12; // Bytes consumed so far
      } else {
        // Format 0, 4, 6 Header: format(u16) + length(u16) + language(u16) + ...
        length = readT(u16, bytes, &sub_cursor);
        language = readT(u16, bytes, &sub_cursor);
        header_size = 6; // Bytes consumed so far
      }

      // Calculate size of the remaining body (the "glyph_index" payload)
      // Safety check to ensure the font file isn't truncated
      if (rec.offset + length > bytes.len) return error.EndOfStream;

      const payload_size = length - header_size;

      const raw_payload = bytes[sub_cursor..][0..payload_size];

      tables[i] = CmapSubtable{
        .format = format,
        .length = @truncate(length),
        .language = language,
        .glyph_index = raw_payload,
      };
    }

    return CmapTable{
      .header = header,
      .records = records,
      .tables = tables,
    };
  }

  pub fn deinit(self: Self, allocator: std.mem.Allocator) void {
    allocator.free(self.records);
    allocator.free(self.tables);
  }

  pub fn pretty_print(self: Self) void {
    const print = std.debug.print;

    print("Cmap Table:\n", .{});
    print("  Header:\n", .{});
    print("    Version: {d}\n", .{self.header.version});
    print("    Num Subtables: {d}\n", .{self.header.numberSubtables});

    print("  Records:\n", .{});
    for (self.records, 0..) |rec, i| {
      print("    [{d}] Platform: {s} ({d})\n", .{
        i, @tagName(rec.platformID), @intFromEnum(rec.platformID)
      });
      print("        Encoding ID: {d}\n", .{rec.platformSpecificID});
      print("        Offset: 0x{x}\n", .{rec.offset});
    }

    print("  Subtables:\n", .{});
    for (self.tables, 0..) |table, i| {
      print("    [{d}] Format {d}\n", .{ i, table.format });
      print("        Length: {d}\n", .{table.length});
      print("        Language: {d}\n", .{table.language});

      // Print a small preview of the raw bytes
      print("        Raw Data: {d} bytes [ ", .{table.glyph_index.len});
      const preview_len = @min(table.glyph_index.len, 8);
      for (table.glyph_index[0..preview_len]) |b| {
          print("{x:0>2} ", .{b});
      }
      if (table.glyph_index.len > preview_len) {
          print("... ", .{});
      }
      print("]\n", .{});
    }
}
};

const Point = struct {
  x: i16,
  y: i16,
  on_curve: bool,
};

pub const GlyphTable = struct {
  number_of_contours: i16,
  x_min: i16,
  y_min: i16,
  x_max: i16,
  y_max: i16,

  points: []Point,
  end_pts_of_contours: []u16, // needed for contour iterator
  instructions: []const u8,

  /// Initialize a simple glyph from bytes
  pub fn init(bytes: []const u8, allocator: *std.mem.Allocator) !GlyphTable {
    var off: usize = 0;

    const number_of_contours = @as(usize, @intCast(readT(i16, bytes, &off)));
    const x_min = readT(i16, bytes, &off);
    const y_min = readT(i16, bytes, &off);
    const x_max = readT(i16, bytes, &off);
    const y_max = readT(i16, bytes, &off);

    if (number_of_contours <= 0) return error.UnsupportedGlyph;

    // --- Read endPtsOfContours ---
    var end_pts_of_contours: []u16 = try allocator.alloc(u16, number_of_contours);
    for (0..number_of_contours) |i| {
      end_pts_of_contours[i] = readT(u16, bytes, &off);
    }

    const total_points = @as(usize, @intCast(end_pts_of_contours[number_of_contours - 1])) + 1;

    // --- Read instructions ---
    const instruction_length: usize = @as(usize, @intCast(readT(u16, bytes, &off)));
    const instructions = bytes[off..off + instruction_length];
    off += instruction_length;

    // --- Read flags ---
    var flags: []u8 = try allocator.alloc(u8, total_points);
    {
      var point_index: usize = 0;
      while (point_index < total_points) {
        const flag: u8 = bytes[off];
        off += 1;
        flags[point_index] = flag;
        point_index += 1;

        if ((flag & 0x08) != 0) { // repeat flag
          const repeat_count: u8 = bytes[off];
          off += 1;
          for (0..repeat_count) |_| {
            if (point_index >= total_points) break;
            flags[point_index] = flag;
            point_index += 1;
          }
        }
      }
    }

    // --- Read xCoordinates ---
    var x_coords: []i16 = try allocator.alloc(i16, total_points);
    var x: i16 = 0;
    for (0..total_points) |point_index| {
      const flag = flags[point_index];
      if ((flag & 0x02) != 0) { // 1-byte x
        const dx: u8 = bytes[off];
        off += 1;
        x += if ((flag & 0x10) != 0) @as(i16, @intCast(dx)) else -@as(i16, @intCast(dx));
      } else if ((flag & 0x10) == 0) { // 2-byte x
        x += readT(i16, bytes, &off);
      }
      x_coords[point_index] = x;
    }

    // --- Read yCoordinates ---
    var y_coords: []i16 = try allocator.alloc(i16, total_points);
    var y: i16 = 0;
    for (0..total_points) |point_index| {
      const flag = flags[point_index];
      if ((flag & 0x04) != 0) { // 1-byte y
        const dy: u8 = bytes[off];
        off += 1;
        y += if ((flag & 0x20) != 0) @as(i16, @intCast(dy)) else -@as(i16, @intCast(dy));
      } else if ((flag & 0x20) == 0) { // 2-byte y
        y += readT(i16, bytes, &off);
      }
      y_coords[point_index] = y;
    }

    // --- Build points array ---
    var points: []Point = try allocator.alloc(Point, total_points);
    for (0..total_points) |point_index| {
      points[point_index] = Point{
        .x = x_coords[point_index],
        .y = y_coords[point_index],
        .on_curve = (flags[point_index] & 0x01) != 0,
      };
    }

    // Free temporary arrays
    allocator.free(flags);
    allocator.free(x_coords);
    allocator.free(y_coords);

    return GlyphTable{
      .number_of_contours = number_of_contours,
      .x_min = x_min,
      .y_min = y_min,
      .x_max = x_max,
      .y_max = y_max,
      .points = points,
      .end_pts_of_contours = end_pts_of_contours,
      .instructions = instructions,
    };
  }

  /// Contour iterator
  pub fn contours(self: *GlyphTable) ContourIterator {
    return ContourIterator{
      .glyph = self,
      .index = 0,
      .start_point = 0,
    };
  }

  pub const ContourIterator = struct {
    glyph: *GlyphTable,
    index: usize,
    start_point: usize,

    pub fn next(self: *ContourIterator) ?[]Point {
      if (self.index >= @as(usize, @intCast(self.glyph.number_of_contours))) return null;

      const end_point = @as(usize, @intCast(self.glyph.end_pts_of_contours[self.index]));
      const contour_points = self.glyph.points[self.start_point .. end_point + 1];

      self.index += 1;
      self.start_point = end_point + 1;

      return contour_points;
    }
  };
};

pub const TtfFont = struct {
  ttfcontent: []const u8,
  offset_table: OffsetTable,
  table_records: []TableRecord,

  pub fn load(allocator: std.mem.Allocator, ttfcontent: []const u8) !TtfFont {
    var offset: usize = 0;
    const offset_table = OffsetTable.fromBytes(ttfcontent, &offset);
    const table_records = try allocator.alloc(TableRecord, offset_table.num_tables);
    std.log.debug("offset {} offset_table {any}", .{ offset, offset_table });
    for (table_records) |*tr| {
      tr.* = TableRecord.fromBytes(ttfcontent, &offset);
      std.log.debug("offset {}", .{ offset });
    }

    return TtfFont{
      .ttfcontent = ttfcontent,
      .offset_table = offset_table,
      .table_records = table_records,
    };
  }

  const Self = @This();
  pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
    allocator.free(self.table_records);
  }

  fn getTableRecords(self: Self, tag: []const u8) ?TableRecord {
    for (self.table_records) |table_record| {
      if (std.mem.eql(u8, &table_record.tag, tag)) {
        return table_record;
      }
    }
    return null;
  }

  pub fn getHead(self: Self) !HeadTable {
    const head_record = self.getTableRecords("head");
    if (head_record) |hr| {
      const head_table = try HeadTable.init(self.ttfcontent[hr.offset..hr.offset + hr.length]);
      return head_table;
    }
    return error.NoHeadTable;
  }

  pub fn getCmap(self: Self, allocator: std.mem.Allocator) !CmapTable {
    const cmap_record = self.getTableRecords("cmap");
    if (cmap_record) |cr| {
      const cmap_table = try CmapTable.init(allocator, self.ttfcontent[cr.offset..cr.offset + cr.length]);
      return cmap_table;
    }
    return error.NoHeadTable;
  }
};
