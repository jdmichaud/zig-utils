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

  fn fromBytes(bytes: []const u8) OffsetTable {
    var off: usize = 0;
    return .{
      .scaler_type = readT(u32, bytes, &off),
      .num_tables = readT(u16, bytes, &off),
      .search_range = readT(u16, bytes, &off),
      .entry_selector = readT(u16, bytes, &off),
      .range_shift = readT(u16, bytes, &off),
    };
  }

};

const TableRecord = struct {
  tag: [4]u8,
  checksum: u32,
  offset: u32,
  length: u32,

  fn fromBytes(bytes: []const u8) TableRecord {
    var tag: [4]u8 = undefined;
    @memcpy(&tag, bytes[0..4]);

    var off: usize = 4;
    return .{
      .tag = tag,
      .checksum = readT(u32, bytes, &off),
      .offset = readT(u32, bytes, &off),
      .length = readT(u32, bytes, &off),
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
    const offset_table = OffsetTable.fromBytes(ttfcontent);
    const table_records = try allocator.alloc(TableRecord, offset_table.num_tables);

    var offset: usize = 12; // size of TableRecord struct
    for (0..offset_table.num_tables) |i| {
      table_records[i] = TableRecord.fromBytes(ttfcontent[offset..]);
      offset += 16;
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
};
