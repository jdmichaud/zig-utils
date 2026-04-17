const std = @import("std");

/// Create a type where each field of T becomes optional with null default.
/// For non-struct types, returns T unchanged.
pub fn Partial(comptime T: type) type {
  return switch (@typeInfo(T)) {
    .@"struct" => |info| blk: {
      var fields: [info.fields.len]std.builtin.Type.StructField = undefined;
      for (info.fields, 0..) |field, i| {
        fields[i] = .{
          .name = field.name,
          .type = ?field.type,
          .default_value_ptr = @ptrCast(&@as(?field.type, null)),
          .is_comptime = false,
          .alignment = 0,
        };
      }
      break :blk @Type(.{
        .@"struct" = .{
          .layout = .auto,
          .fields = &fields,
          .decls = &.{},
          .is_tuple = false,
        },
      });
    },
    else => T,
  };
}

/// Check whether `value` matches `pattern`.
/// For struct types: null fields in pattern are wildcards.
/// For non-struct types: exact equality.
pub fn partialMatches(comptime T: type, pattern: Partial(T), value: T) bool {
  switch (@typeInfo(T)) {
    .@"struct" => |info| {
      inline for (info.fields) |field| {
        if (@field(pattern, field.name)) |val| {
          if (val != @field(value, field.name)) return false;
        }
      }
      return true;
    },
    else => return std.meta.eql(pattern, value),
  }
}
