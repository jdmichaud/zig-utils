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

pub const MaxpTable = struct {
  /// 0x00010000 (1.0) or 0x00005000 (0.5)
  version: Fixed,
  /// The number of glyphs in the font.
  num_glyphs: u16,

  // The following are only present if version is 1.0
  max_points: u16 = 0,
  max_contours: u16 = 0,
  max_composite_points: u16 = 0,
  max_composite_contours: u16 = 0,
  max_zones: u16 = 0,
  max_twilight_points: u16 = 0,
  max_storage: u16 = 0,
  max_function_defs: u16 = 0,
  max_instruction_defs: u16 = 0,
  max_stack_elements: u16 = 0,
  max_size_of_instructions: u16 = 0,
  max_component_elements: u16 = 0,
  max_component_depth: u16 = 0,

  const Self = @This();

  pub fn pretty_print(self: Self) void {
    const print = std.debug.print;
    print("Maxp Table:\n", .{});
    print("  version: {}\n", .{self.version});
    print("  num_glyphs: {}\n", .{self.num_glyphs});

    // Check if version is 1.0 (0x00010000)
    if (self.version.raw == 0x00010000) {
      print("  max_points: {}\n", .{self.max_points});
      print("  max_contours: {}\n", .{self.max_contours});
      print("  max_composite_points: {}\n", .{self.max_composite_points});
      print("  max_composite_contours: {}\n", .{self.max_composite_contours});
      print("  max_zones: {}\n", .{self.max_zones});
      print("  max_twilight_points: {}\n", .{self.max_twilight_points});
      print("  max_storage: {}\n", .{self.max_storage});
      print("  max_function_defs: {}\n", .{self.max_function_defs});
      print("  max_instruction_defs: {}\n", .{self.max_instruction_defs});
      print("  max_stack_elements: {}\n", .{self.max_stack_elements});
      print("  max_size_of_instructions: {}\n", .{self.max_size_of_instructions});
      print("  max_component_elements: {}\n", .{self.max_component_elements});
      print("  max_component_depth: {}\n", .{self.max_component_depth});
    }
  }

  pub fn init(bytes: []const u8) !MaxpTable {
    std.log.debug("bytes.len {}", .{ bytes.len });
    // Minimum size for version 0.5 is 6 bytes
    if (bytes.len < 6) return error.UnexpectedEndOfTable;

    var off: usize = 0;
    const version_val = readT(i32, bytes, &off);
    const version = Fixed{ .raw = version_val };
    const num_glyphs = readT(u16, bytes, &off);

    // Basic struct with defaults
    var table = MaxpTable{
      .version = version,
      .num_glyphs = num_glyphs,
    };

    // If version is 1.0 (0x00010000), parse the extra fields
    if (version_val == 0x00010000) {
      if (bytes.len < 32) return error.UnexpectedEndOfTable;

      table.max_points = readT(u16, bytes, &off);
      table.max_contours = readT(u16, bytes, &off);
      table.max_composite_points = readT(u16, bytes, &off);
      table.max_composite_contours = readT(u16, bytes, &off);
      table.max_zones = readT(u16, bytes, &off);
      table.max_twilight_points = readT(u16, bytes, &off);
      table.max_storage = readT(u16, bytes, &off);
      table.max_function_defs = readT(u16, bytes, &off);
      table.max_instruction_defs = readT(u16, bytes, &off);
      table.max_stack_elements = readT(u16, bytes, &off);
      table.max_size_of_instructions = readT(u16, bytes, &off);
      table.max_component_elements = readT(u16, bytes, &off);
      table.max_component_depth = readT(u16, bytes, &off);
    }

    return table;
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

  const Self = @This();
  pub fn pretty_print(self: Self) void {
    const print = std.debug.print;

    print("Head Table:\n", .{});
    print("  version {} \n", .{ self.version });
    print("  font_revision {} \n", .{ self.font_revision });
    // print("  flags {} \n", .{ self.flags });
    print("  units_per_em {} \n", .{ self.units_per_em });

    print("  x_min {} \n", .{ self.x_min });
    print("  y_min {} \n", .{ self.y_min });
    print("  x_max {} \n", .{ self.x_max });
    print("  y_max {} \n", .{ self.y_max });

    print("  index_to_loc_format {} \n", .{ self.index_to_loc_format });
    print("  glyph_data_format {} \n", .{ self.glyph_data_format });
  }

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

pub const HheaTable = struct {
  /// 16.16 fixed-point number
  /// Major version number of the horizontal header table, set to 1.0 (0x00010000)
  version: Fixed,

  /// Typographic ascent (distance from baseline of highest ascender)
  ascender: i16,
  /// Typographic descent (distance from baseline of lowest descender)
  descender: i16,
  /// Typographic line gap. Negative LineGap values are treated as zero in some legacy platform implementations.
  line_gap: i16,

  /// Maximum advance width value in 'hmtx' table.
  advance_width_max: u16,
  /// Minimum Left Side Bearing value in 'hmtx' table.
  min_left_side_bearing: i16,
  /// Minimum Right Side Bearing value; calculated as (advanceWidth - lsb - (xMax - xMin)).
  min_right_side_bearing: i16,
  /// Max(lsb + (xMax - xMin)).
  x_max_extent: i16,

  /// Used to calculate the slope of the cursor (rise/run); 1 for vertical.
  caret_slope_rise: i16,
  /// 0 for vertical.
  caret_slope_run: i16,
  /// The amount by which a slanted highlight on a glyph needs to be shifted to produce the best appearance. Set to 0 for non-slanted fonts.
  caret_offset: i16,

  /// 0 for current format.
  metric_data_format: i16,
  /// Number of hMetric entries in 'hmtx' table
  number_of_h_metrics: u16,

  const Self = @This();

  pub fn pretty_print(self: Self) void {
    const print = std.debug.print;

    print("Hhea Table:\n", .{});
    print("  version {} \n", .{self.version});
    print("  ascender {} \n", .{self.ascender});
    print("  descender {} \n", .{self.descender});
    print("  line_gap {} \n", .{self.line_gap});

    print("  advance_width_max {} \n", .{self.advance_width_max});
    print("  min_left_side_bearing {} \n", .{self.min_left_side_bearing});
    print("  min_right_side_bearing {} \n", .{self.min_right_side_bearing});
    print("  x_max_extent {} \n", .{self.x_max_extent});

    print("  caret_slope_rise {} \n", .{self.caret_slope_rise});
    print("  caret_slope_run {} \n", .{self.caret_slope_run});
    print("  caret_offset {} \n", .{self.caret_offset});

    print("  metric_data_format {} \n", .{self.metric_data_format});
    print("  number_of_h_metrics {} \n", .{self.number_of_h_metrics});
  }

  // Size of the hhea table (36 bytes)
  pub const byte_len: usize = 36;

  /// Parse an `hhea` table from a big-endian byte slice
  pub fn init(bytes: []const u8) !HheaTable {
    if (bytes.len < byte_len) {
        return error.UnexpectedEndOfTable;
    }

    var off: usize = 0;

    const version = Fixed{ .raw = readT(i32, bytes, &off) };
    const ascender = readT(i16, bytes, &off);
    const descender = readT(i16, bytes, &off);
    const line_gap = readT(i16, bytes, &off);

    const advance_width_max = readT(u16, bytes, &off);
    const min_left_side_bearing = readT(i16, bytes, &off);
    const min_right_side_bearing = readT(i16, bytes, &off);
    const x_max_extent = readT(i16, bytes, &off);

    const caret_slope_rise = readT(i16, bytes, &off);
    const caret_slope_run = readT(i16, bytes, &off);
    const caret_offset = readT(i16, bytes, &off);

    // Skip 4 reserved int16 fields (8 bytes total)
    _ = readT(i16, bytes, &off);
    _ = readT(i16, bytes, &off);
    _ = readT(i16, bytes, &off);
    _ = readT(i16, bytes, &off);

    const metric_data_format = readT(i16, bytes, &off);
    const number_of_h_metrics = readT(u16, bytes, &off);

    return HheaTable{
      .version = version,
      .ascender = ascender,
      .descender = descender,
      .line_gap = line_gap,
      .advance_width_max = advance_width_max,
      .min_left_side_bearing = min_left_side_bearing,
      .min_right_side_bearing = min_right_side_bearing,
      .x_max_extent = x_max_extent,
      .caret_slope_rise = caret_slope_rise,
      .caret_slope_run = caret_slope_run,
      .caret_offset = caret_offset,
      .metric_data_format = metric_data_format,
      .number_of_h_metrics = number_of_h_metrics,
    };
  }
};

// Origin (x=0)                                        Next Origin
//       |                                          (Where the next char starts)
//       v                                                   v
//       |------------------ Advance Width ------------------|
//       |                                                   |
//       |   LSB    |        Glyph Bounding Box      |  RSB  |
//       | <------> | <----------------------------> | <---> |
//       |          |________________________________|       |
//       |          |                                |       |
//       |          |               /\               |       |
//       |          |              /  \              |       |
//       |          |             /    \             |       |
//       |          |            /      \            |       |
//       |          |           /________\           |       |
//       |          |          /          \          |       |
//       |          |         /            \         |       |
//       |__________|________/______________\________|_______|
//                  ^
//                  xMin (from 'head' or glyph data)
// This pseudo code demonstrate how characters are spaced relative to each other
// ```
// pen_x = 0
//
// for each glyph:
//     draw glyph at:
//         x = pen_x + leftSideBearing
//     pen_x += advanceWidth
//     pen_x += kerning(previous_glyph, current_glyph)
// ```
// `advanceWidth` and `leftSideBearing` is to be found in `hmtx` table
// `kerning` is to be found in `kern` or `GPOS` table
// `hhea` table provides information on how to parse the `hmtx` table
pub const LongHorMetric = struct {
  /// The amount of displacement in X after a character
  advance_width: u16,
  /// left side bearing is the distance from the left bound to the left most contour
  lsb: i16,
};

pub const HmtxTable = struct {
  number_of_h_metrics: u16,
  num_glyphs: u16,

  /// Slice containing the hMetrics array (numberOfHMetrics * 4 bytes)
  h_metrics_data: []const u8,
  /// Slice containing the leftSideBearing array (remaining glyphs * 2 bytes)
  lsb_data: []const u8,

  const Self = @This();

  pub fn pretty_print(self: Self) void {
    const print = std.debug.print;
    print("Hmtx Table:\n", .{});
    print("  number_of_h_metrics {} \n", .{self.number_of_h_metrics});
    print("  num_glyphs {} \n", .{self.num_glyphs});
    print("  data size: {} bytes\n", .{self.h_metrics_data.len + self.lsb_data.len});
  }

  /// Parse an `hmtx` table.
  /// Note: This table requires `number_of_h_metrics` from the `hhea` table
  /// and `num_glyphs` from the `maxp` table to be parsed correctly.
  pub fn init(bytes: []const u8, number_of_h_metrics: u16, num_glyphs: u16) !HmtxTable {
    // Validation: number_of_h_metrics should not exceed num_glyphs
    if (number_of_h_metrics > num_glyphs) {
      // While technically the spec says this shouldn't happen, some fonts might be malformed.
      // A robust parser might clamp this or error. Here we return error.
      return error.InvalidHmtxCounts;
    }

    var off: usize = 0;

    // 1. hMetrics array: numberOfHMetrics * 4 bytes
    // Each entry is { advanceWidth: u16, lsb: i16 }
    const h_metrics_len = @as(usize, number_of_h_metrics) * 4;
    if (off + h_metrics_len > bytes.len) return error.UnexpectedEndOfTable;

    const h_metrics_slice = bytes[off .. off + h_metrics_len];
    off += h_metrics_len;

    // 2. leftSideBearing array: (numGlyphs - numberOfHMetrics) * 2 bytes
    // Each entry is { lsb: i16 }
    // The advanceWidth for these glyphs is the same as the last entry in hMetrics.
    const num_lsbs = num_glyphs - number_of_h_metrics;
    const lsb_len = @as(usize, num_lsbs) * 2;
    if (off + lsb_len > bytes.len) return error.UnexpectedEndOfTable;

    const lsb_slice = bytes[off .. off + lsb_len];

    return HmtxTable{
      .number_of_h_metrics = number_of_h_metrics,
      .num_glyphs = num_glyphs,
      .h_metrics_data = h_metrics_slice,
      .lsb_data = lsb_slice,
    };
  }

  /// Retrieve the horizontal metrics for a specific glyph index.
  pub fn getMetric(self: Self, glyph_index: u16) !LongHorMetric {
    if (glyph_index >= self.num_glyphs) {
      return error.GlyphIndexOutOfBounds;
    }

    if (glyph_index < self.number_of_h_metrics) {
      // Case 1: Glyph is inside the hMetrics array.
      // We read both advance_width and lsb from the first array.
      var off = @as(usize, glyph_index) * 4;
      return LongHorMetric{
        .advance_width = readT(u16, self.h_metrics_data, &off),
        .lsb = readT(i16, self.h_metrics_data, &off),
      };
    } else {
      // Case 2: Glyph is inside the leftSideBearing array (monospaced remainder).
      // Its advance_width is the same as the LAST entry in the hMetrics array.

      // Get last width
      const last_metric_idx = self.number_of_h_metrics - 1;
      var width_off = @as(usize, last_metric_idx) * 4;
      const common_width = readT(u16, self.h_metrics_data, &width_off);

      // Get specific LSB
      const lsb_index = glyph_index - self.number_of_h_metrics;
      var lsb_off = @as(usize, lsb_index) * 2;
      const lsb = readT(i16, self.lsb_data, &lsb_off);

      return LongHorMetric{
        .advance_width = common_width,
        .lsb = lsb,
      };
    }
  }
};

// cmap maps character codes (e.g. Unicode) to glyph IDs (GIDs). The GID indexes
// into the `offsets` field used by loca. loca stores offsets (or indices) into
// the glyf table relative to the beginning of the glyf table for each GID. Each
// loca entry gives the start (and next gives the end) of a glyph. glyf contains
// the actual glyph data (outlines, components, hints).
// Flow: char code → cmap → GID → loca → glyf.
// In other words: glyf[loca.offsets[cmap[code_point]]]

const CmapHeader = struct {
  /// Version number (Set to zero)
  version: u16,
  /// Number of encoding subtables
  number_sub_tables: u16,
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

/// The specific implementation for Format 4 (Standard BMP)
const CmapFormat4 = struct {
  seg_count: u16,
  // We store pre-calculated SLICES.
  // This means we don't have to do "offset = 14 + i * 2" math every time.
  end_counts: []const u8,
  start_counts: []const u8,
  id_deltas: []const u8,
  id_range_offsets: []const u8,
  glyph_id_array: []const u8, // The rest of the data

  // We need the full subtable data for idRangeOffset relative math
  raw_subtable: []const u8,
  idRangeOffset_start_pos: usize,

  pub fn init(bytes: []const u8) !CmapFormat4 {
    if (bytes.len < 14) return error.InvalidCmapTable;

    const seg_countX2 = std.mem.readInt(u16, bytes[6..8], .big);
    const seg_count = seg_countX2 / 2;

    // Calculate offsets ONCE
    const end_cnt_off = 14;
    const reserved_pad = end_cnt_off + seg_countX2;
    const start_cnt_off = reserved_pad + 2;
    const id_delta_off = start_cnt_off + seg_countX2;
    const id_range_off = id_delta_off + seg_countX2;
    const glyph_id_arr_off = id_range_off + seg_countX2;

    if (glyph_id_arr_off > bytes.len) return error.InvalidCmapTable;

    return CmapFormat4{
      .seg_count = seg_count,
      .end_counts = bytes[end_cnt_off..][0..seg_countX2],
      .start_counts = bytes[start_cnt_off..][0..seg_countX2],
      .id_deltas = bytes[id_delta_off..][0..seg_countX2],
      .id_range_offsets = bytes[id_range_off..][0..seg_countX2],
      .glyph_id_array = bytes[glyph_id_arr_off..], // The rest
      .raw_subtable = bytes,
      .idRangeOffset_start_pos = id_range_off,
    };
  }

  pub fn get(self: CmapFormat4, cp: u32) u16 {
    if (cp > 65535) return 0;
    const cp_u16 = @as(u16, @truncate(cp));

    // 1. Find Segment (Linear search is fast enough for small seg_counts,
    // but Binary Search is better for full fonts. Kept simple here.)
    var i: u16 = 0;
    while (i < self.seg_count) : (i += 1) {
      // We read directly from our pre-sliced array
      const end = std.mem.readInt(u16, self.end_counts[i * 2..][0..2], .big);
      if (end >= cp_u16) break;
    }
    if (i == self.seg_count) return 0;

    const start = std.mem.readInt(u16, self.start_counts[i * 2..][0..2], .big);
    if (start > cp_u16) return 0;

    const id_delta = std.mem.readInt(u16, self.id_deltas[i * 2..][0..2], .big);
    const id_range_offset = std.mem.readInt(u16, self.id_range_offsets[i * 2..][0..2], .big);

    if (id_range_offset == 0) {
      return @as(u16, @truncate(cp_u16 +% id_delta));
    }

    // Pointer math logic
    // Location of the specific rangeOffset we just read
    const current_iro_loc = self.idRangeOffset_start_pos + (i * 2);

    // Offset relative to that location
    const offset = current_iro_loc + id_range_offset + (cp_u16 - start) * 2;

    if (offset + 2 > self.raw_subtable.len) return 0;

    const glyph = std.mem.readInt(u16, self.raw_subtable[offset..][0..2], .big);
    if (glyph == 0) return 0;

    return @as(u16, @truncate(glyph +% id_delta));
  }

  /// Calculates the total number of characters mapped by this subtable
  pub fn getMappedCharCount(self: CmapFormat4) u32 {
    var total: u32 = 0;
    var i: u16 = 0;
    while (i < self.seg_count) : (i += 1) {
      const start = std.mem.readInt(u16, self.start_counts[i * 2..][0..2], .big);
      const end = std.mem.readInt(u16, self.end_counts[i * 2..][0..2], .big);

      // Add the range size (inclusive)
      if (end >= start) {
        total += (end - start) + 1;
      }
    }
    return total;
  }
};

/// The specific implementation for Format 12 (Modern)
const CmapFormat12 = struct {
  num_groups: u32,
  groups_data: []const u8, // Slice starting at the first group

  pub fn init(bytes: []const u8) !CmapFormat12 {
    if (bytes.len < 16) return error.InvalidCmapTable;

    const num_groups = std.mem.readInt(u32, bytes[12..16], .big);

    // Verify we have enough bytes for all groups (12 bytes per group)
    if (16 + num_groups * 12 > bytes.len) return error.InvalidCmapTable;

    return CmapFormat12{
      .num_groups = num_groups,
      .groups_data = bytes[16..],
    };
  }

  pub fn get(self: CmapFormat12, cp: u32) u16 {
    // Binary search is ideal here, but linear for clarity:
    var i: u32 = 0;
    var offset: usize = 0;
    while (i < self.num_groups) : (i += 1) {
      const start = std.mem.readInt(u32, self.groups_data[offset..][0..4], .big);
      const end = std.mem.readInt(u32, self.groups_data[offset + 4..][0..4], .big);

      if (cp >= start and cp <= end) {
        const start_id = std.mem.readInt(u32, self.groups_data[offset + 8..][0..4], .big);
        // Although format 12 stored 32 bits, glyph index is a 16 bits. The maxp
        // table `num_glpyh` is a u16. That is why TTF font do not support more
        // than 65K characters.
        return @intCast(start_id + (cp - start));
      }
      offset += 12;
    }
    return 0;
  }

  /// Calculates the total number of characters mapped by this subtable
  pub fn getMappedCharCount(self: CmapFormat12) u32 {
    var total: u32 = 0;
    var i: u32 = 0;
    var offset: usize = 0;
    while (i < self.num_groups) : (i += 1) {
      const start = std.mem.readInt(u32, self.groups_data[offset..][0..4], .big);
      const end = std.mem.readInt(u32, self.groups_data[offset + 4..][0..4], .big);

      total += (end - start) + 1;

      offset += 12;
    }
    return total;
  }};

/// The Container that holds "One of the specific implementations"
pub const CmapMapper = union(enum) {
  fmt4: CmapFormat4,
  fmt12: CmapFormat12,
  unknown: void,

  /// This is the "Factory" function
  pub fn init(bytes: []const u8) !CmapMapper {
    if (bytes.len < 2) return error.InvalidCmapTable;
    const format = std.mem.readInt(u16, bytes[0..2], .big);

    switch (format) {
      4 => return CmapMapper{ .fmt4 = try CmapFormat4.init(bytes) },
      12 => return CmapMapper{ .fmt12 = try CmapFormat12.init(bytes) },
      else => return CmapMapper{ .unknown = {} },
    }
  }

  /// The single unified call point
  pub fn getGlyph(self: CmapMapper, code_point: u32) !u16 {
    switch (self) {
      .fmt4 => |impl| return impl.get(code_point),
      .fmt12 => |impl| return impl.get(code_point),
      .unknown => return error.UnknownFormat,
    }
  }

  pub fn getMappedCharCount(self: CmapMapper) !u32 {
    switch (self) {
      .fmt4 => |impl| return impl.getMappedCharCount(),
      .fmt12 => |impl| return impl.getMappedCharCount(),
      .unknown => return error.UnknownFormat,
    }
  }
};

const CmapSubtable = struct {
  format: u16,
  /// length in bytes including this structure
  length: u16,
  /// Language code
  language: u16,
  /// The object that will fetch the glyph ID depending on the format of the subtable
  mapper: CmapMapper,
};

pub const CmapTable = struct {
  header: CmapHeader,
  records: []CmapRecord,
  sub_tables: []CmapSubtable,

  const Self = @This();
  pub fn init(allocator: std.mem.Allocator, bytes: []const u8) !CmapTable {
    var index: usize = 0;

    // 1. Read the Header
    const header = CmapHeader{
      .version = readT(u16, bytes, &index),
      .number_sub_tables = readT(u16, bytes, &index),
    };

    // 2. Allocate and Read Encoding Records
    const records = try allocator.alloc(CmapRecord, header.number_sub_tables);
    errdefer allocator.free(records);
    for (records) |*rec| {
      rec.platformID = @enumFromInt(readT(u16, bytes, &index));
      rec.platformSpecificID = readT(u16, bytes, &index);
      rec.offset = readT(u32, bytes, &index);
    }

    // 3. Read Subtables
    // Note: The records contain offsets relative to the start of the Cmap Table.
    // We must jump to those offsets to read the actual subtable data.
    const sub_tables = try allocator.alloc(CmapSubtable, header.number_sub_tables);
    errdefer allocator.free(sub_tables);

    for (records, 0..) |rec, i| {
      // Create a cursor at the offset for this specific subtable
      var sub_cursor = @as(usize, rec.offset);
      const subtable_bytes = bytes[sub_cursor..];

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

      sub_tables[i] = CmapSubtable{
        .format = format,
        .length = @truncate(length),
        .language = language,
        .mapper = try CmapMapper.init(subtable_bytes),
      };
    }

    return CmapTable{
      .header = header,
      .records = records,
      .sub_tables = sub_tables,
    };
  }

  pub fn deinit(self: Self, allocator: std.mem.Allocator) void {
    allocator.free(self.records);
    allocator.free(self.sub_tables);
  }

  pub fn pretty_print(self: Self) void {
    const print = std.debug.print;

    print("Cmap Table:\n", .{});
    print("  Header:\n", .{});
    print("    Version: {d}\n", .{self.header.version});
    print("    Num Subtables: {d}\n", .{self.header.number_sub_tables});

    print("  Records:\n", .{});
    for (self.records, 0..) |rec, i| {
      print("    [{d}] Platform: {s} ({d})\n", .{
        i, @tagName(rec.platformID), @intFromEnum(rec.platformID)
      });
      print("        Encoding ID: {d}\n", .{rec.platformSpecificID});
      print("        Offset: 0x{x}\n", .{rec.offset});
    }

    print("  Subtables:\n", .{});
    for (self.sub_tables, 0..) |table, i| {
      print("    [{d}] Format {d}\n", .{ i, table.format });
      print("        Length: {d} bytes\n", .{table.length});
      print("        Language: {d}\n", .{table.language});

      // Print a small preview of the raw bytes
      // switch (table.mapper.getMappedCharCount()) {
      //   error => print("        Unknown\n", .{}),
      //   else => |value| print("        Glyphs count: {} glyphs\n", .{ value }),
      // }
      if (table.mapper.getMappedCharCount()) |value| {
        print("        Glyphs count: {} glyphs\n", .{ value });
      } else |_| {
        print("        Unknown format\n", .{});
      }
    }
  }
};

pub const LocaTable = struct {
    /// Offsets of glyphs in the glyphs table in big endian.
    offsets: []const u8,
    is_short: bool,

    const Self = @This();

    pub fn init(index_to_loc_format: i16, bytes: []const u8) !Self {
      // Validate format (Standard says 0 or 1)
      if (index_to_loc_format != 0 and index_to_loc_format != 1) {
        return error.InvalidLocaFormat;
      }

      const is_short = (index_to_loc_format == 0);

      // Optional: Sanity check length alignment
      // Short (u16) needs 2-byte alignment, Long (u32) needs 4-byte
      const alignment = if (is_short) @as(usize, 2) else 4;
      if (bytes.len % alignment != 0) return error.MalformedLocaTable;

      return Self{
        .offsets = bytes,
        .is_short = is_short,
      };
    }

    pub fn get(self: Self, glyph_id: usize) !u32 {
      if (self.is_short) {
        // --- Format 0 (Short) ---
        // Stride is 2 bytes
        const pos = glyph_id * 2;
        if (pos + 2 > self.offsets.len) return error.InvalidGlyphId;

        // Read u16 and multiply by 2
        const val = std.mem.readInt(u16, self.offsets[pos..][0..2], .big);
        return @as(u32, val) * 2;
      } else {
        // --- Format 1 (Long) ---
        // Stride is 4 bytes
        const pos = glyph_id * 4;
        if (pos + 4 > self.offsets.len) return error.InvalidGlyphId;

        // Read u32 directly
        return std.mem.readInt(u32, self.offsets[pos..][0..4], .big);
      }
    }

    pub fn count(self: Self) usize {
      if (self.is_short) {
        return self.offsets.len / 2;
      } else {
        return self.offsets.len / 4;
      }
  }
};

const Point = struct {
  x: i16,
  y: i16,
  on_curve: bool,
};

pub const GlyphTable = struct {
  bytes: []const u8,

  const Self = @This();
  pub fn init(bytes: []const u8) Self {
    return Self{
      .bytes = bytes,
    };
  }

  pub fn getGlyph(self: Self, allocator: std.mem.Allocator, code_point: u32,
    cmap_sub_table: CmapSubtable, loca_table: LocaTable) !Glyph {
    const offset = try loca_table.get(try cmap_sub_table.mapper.getGlyph(code_point));
    if (offset >= self.bytes.len) {
      return error.InvalidGlyphOffset;
    }

    return Glyph.init(allocator, self.bytes[offset..]);
  }
};

pub const Glyph = struct {
  number_of_contours: i16,
  x_min: i16,
  y_min: i16,
  x_max: i16,
  y_max: i16,

  points: []Point,
  end_pts_of_contours: []u16, // needed for contour iterator
  instructions: []const u8,

  const Self = @This();
  /// Initialize a simple glyph from bytes
  pub fn init(allocator: std.mem.Allocator, bytes: []const u8) !Self {
    var off: usize = 0;

    const number_of_contours = readT(i16, bytes, &off);

    const number_of_contours_usize = @as(usize, @intCast(number_of_contours));
    const x_min = readT(i16, bytes, &off);
    const y_min = readT(i16, bytes, &off);
    const x_max = readT(i16, bytes, &off);
    const y_max = readT(i16, bytes, &off);

    if (number_of_contours_usize <= 0) return error.UnsupportedGlyph;

    // --- Read endPtsOfContours ---
    var end_pts_of_contours: []u16 = try allocator.alloc(u16, number_of_contours_usize);
    errdefer allocator.free(end_pts_of_contours);
    for (0..number_of_contours_usize) |i| {
      end_pts_of_contours[i] = readT(u16, bytes, &off);
    }

    const total_points = @as(usize, @intCast(end_pts_of_contours[number_of_contours_usize - 1])) + 1;

    // --- Read instructions ---
    const instruction_length: usize = @as(usize, @intCast(readT(u16, bytes, &off)));
    const instructions = bytes[off..off + instruction_length];
    off += instruction_length;

    // --- Read flags ---
    var flags: []u8 = try allocator.alloc(u8, total_points);
    defer allocator.free(flags);
    {
      var point_index: usize = 0;
      while (point_index < total_points) {
        const flag = readT(u8, bytes, &off);
        flags[point_index] = flag;
        point_index += 1;

        if ((flag & 0x08) != 0) { // repeat flag
          const repeat_count = readT(u8, bytes, &off);
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
    defer allocator.free(x_coords);
    var x: i16 = 0;
    for (0..total_points) |point_index| {
      const flag = flags[point_index];
      if ((flag & 0x02) != 0) { // 1-byte x
        const dx = readT(u8, bytes, &off);
        x += if ((flag & 0x10) != 0) @as(i16, @intCast(dx)) else -@as(i16, @intCast(dx));
      } else if ((flag & 0x10) == 0) { // 2-byte x
        x += readT(i16, bytes, &off);
      }
      x_coords[point_index] = x;
    }

    // --- Read yCoordinates ---
    var y_coords: []i16 = try allocator.alloc(i16, total_points);
    defer allocator.free(y_coords);
    var y: i16 = 0;
    for (0..total_points) |point_index| {
      const flag = flags[point_index];
      if ((flag & 0x04) != 0) { // 1-byte y
        const dy = readT(u8, bytes, &off);
        y += if ((flag & 0x20) != 0) @as(i16, @intCast(dy)) else -@as(i16, @intCast(dy));
      } else if ((flag & 0x20) == 0) { // 2-byte y
        y += readT(i16, bytes, &off);
      }
      y_coords[point_index] = y;
    }

    // --- Build points array ---
    var points: []Point = try allocator.alloc(Point, total_points);
    errdefer allocator.free(points);
    for (0..total_points) |point_index| {
      points[point_index] = Point{
        .x = x_coords[point_index],
        .y = y_coords[point_index],
        .on_curve = (flags[point_index] & 0x01) != 0,
      };
    }

    return Glyph {
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

  pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
    allocator.free(self.end_pts_of_contours);
    allocator.free(self.points);
  }

  /// Contour iterator
  pub fn contours(self: *Glyph) ContourIterator {
    return ContourIterator{
      .glyph = self,
      .index = 0,
      .start_point = 0,
    };
  }

  pub const ContourIterator = struct {
    glyph: *Glyph,
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

  pub const Segment = union(enum) {
    move_to: Point,
    line_to: Point,
    quad_to: struct { cx: i16, cy: i16, x: i16, y: i16 },
  };

  /// Returns an iterator that yields drawing segments (Lines, Bezier curves)
  pub fn curves(self: *Glyph) CurveIterator {
    return CurveIterator{ .contour_iter = self.contours() };
  }

  pub const CurveIterator = struct {
    contour_iter: ContourIterator,
    current_contour: ?[]Point = null,
    start_index: usize = 0,
    processed_count: usize = 0,
    state: enum { Init, Processing } = .Init,

    pub fn next(self: *CurveIterator) ?Segment {
      while (true) {
        // Load next contour if needed
        if (self.current_contour == null) {
          self.current_contour = self.contour_iter.next() orelse return null;
          self.state = .Init;
        }

        const points = self.current_contour.?;
        const len = points.len;
        if (len == 0) { // Safety check
          self.current_contour = null;
          continue;
        }

        if (self.state == .Init) {
          // Determine the visual start point.
          // In TTF, the first point in the array isn't always the start if it is off-curve.
          var first_on_curve: ?usize = null;
          for (points, 0..) |p, i| {
            if (p.on_curve) {
              first_on_curve = i;
              break;
            }
          }

          // If no on-curve point exists, start loop effectively at len-1 so logic generates midpoints
          self.start_index = first_on_curve orelse (len - 1);
          self.processed_count = 0;
          self.state = .Processing;

          if (first_on_curve) |idx| {
            return Segment{ .move_to = points[idx] };
          } else {
            // Edge case: All points are off-curve. Start at midpoint of first and last.
            const p0 = points[0];
            const pl = points[len - 1];
            const mx = @divTrunc(p0.x + pl.x, 2);
            const my = @divTrunc(p0.y + pl.y, 2);
            return Segment{ .move_to = Point{ .x = mx, .y = my, .on_curve = true } };
          }
        }

        // Processing loop
        if (self.processed_count >= len) {
          self.current_contour = null;
          continue;
        }

        const curr_idx = (self.start_index + 1 + self.processed_count) % len;
        const p = points[curr_idx];

        if (p.on_curve) {
          self.processed_count += 1;
          return Segment{ .line_to = p };
        } else {
          // p is off-curve (control point). Check next point.
          const next_idx = (curr_idx + 1) % len;
          const next_p = points[next_idx];

          if (next_p.on_curve) {
            self.processed_count += 2; // Consumed control and endpoint
            return Segment{ .quad_to = .{ .cx = p.x, .cy = p.y, .x = next_p.x, .y = next_p.y } };
          } else {
            // Next is also off-curve: Implied on-curve point at midpoint
            const mx = @divTrunc(p.x + next_p.x, 2);
            const my = @divTrunc(p.y + next_p.y, 2);
            self.processed_count += 1; // Consumed p, but next_p remains for next segment
            return Segment{ .quad_to = .{ .cx = p.x, .cy = p.y, .x = mx, .y = my } };
          }
        }
      }
    }
  };
};

pub const TtfFont = struct {
  ttfcontent: []const u8,
  offset_table: OffsetTable,
  table_records: []TableRecord,

  maxp_table: MaxpTable,
  head_table: HeadTable,
  hhea_table: HheaTable,
  cmap_table: CmapTable,
  loca_table: LocaTable,
  hmtx_table: HmtxTable,
  glyf_table: GlyphTable,

  pub fn load(allocator: std.mem.Allocator, ttfcontent: []const u8) !TtfFont {
    var offset: usize = 0;
    const offset_table = OffsetTable.fromBytes(ttfcontent, &offset);
    const table_records = try allocator.alloc(TableRecord, offset_table.num_tables);

    for (table_records) |*tr| {
      tr.* = TableRecord.fromBytes(ttfcontent, &offset);
    }

    const mr = getTableRecords(table_records, "maxp") orelse return error.NoMaxpTable;
    const maxp_table = try MaxpTable.init(ttfcontent[mr.offset..mr.offset + mr.length]);
    const hr = getTableRecords(table_records, "head") orelse return error.NoHeadTable;
    const head_table = try HeadTable.init(ttfcontent[hr.offset..hr.offset + hr.length]);
    const her = getTableRecords(table_records, "hhea") orelse return error.NoHeeaTable;
    const hhea_table = try HheaTable.init(ttfcontent[her.offset..her.offset + her.length]);
    const cr = getTableRecords(table_records, "cmap") orelse return error.NoCmapTable;
    const cmap_table = try CmapTable.init(allocator, ttfcontent[cr.offset..cr.offset + cr.length]);
    const lr = getTableRecords(table_records, "loca") orelse return error.NoLocaTable;
    const loca_table = try LocaTable.init(head_table.index_to_loc_format, ttfcontent[lr.offset..lr.offset + lr.length]);
    const hmtxr = getTableRecords(table_records, "hmtx") orelse return error.NoGlyphTable;
    const hmtx_table = try HmtxTable.init(ttfcontent[hmtxr.offset..hmtxr.offset + hmtxr.length],
      hhea_table.number_of_h_metrics, maxp_table.num_glyphs);
    const gr = getTableRecords(table_records, "glyf") orelse return error.NoGlyphTable;
    const glyf_table = GlyphTable.init(ttfcontent[gr.offset..gr.offset + gr.length]);

    return TtfFont{
      .ttfcontent = ttfcontent,
      .offset_table = offset_table,
      .table_records = table_records,
      .maxp_table = maxp_table,
      .head_table = head_table,
      .hhea_table = hhea_table,
      .cmap_table = cmap_table,
      .loca_table = loca_table,
      .hmtx_table = hmtx_table,
      .glyf_table = glyf_table,
    };
  }

  const Self = @This();
  pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
    allocator.free(self.table_records);
    self.cmap_table.deinit(allocator);
  }

  fn getTableRecords(table_records: []const TableRecord, tag: []const u8) ?TableRecord {
    for (table_records) |table_record| {
      if (std.mem.eql(u8, &table_record.tag, tag)) {
        return table_record;
      }
    }
    return null;
  }

  pub fn getGlyph(self: Self, allocator: std.mem.Allocator, code_point: u32) !Glyph {
    return self.glyf_table.getGlyph(allocator, code_point, self.cmap_table.sub_tables[0], self.loca_table);
  }

  pub fn getGlyphMetric(self: Self, code_point: u32) !LongHorMetric {
    const glyph_index = try self.cmap_table.sub_tables[0].mapper.getGlyph(code_point);
    return self.hmtx_table.getMetric(glyph_index);
  }
};
