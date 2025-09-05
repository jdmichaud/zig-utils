const std = @import("std");

// Load a file into a buffer using mmap
pub fn load(pathname: []const u8) ![]align(4096) const u8 {
  var file = try std.fs.cwd().openFile(pathname, .{});
  defer file.close();

  const size = try file.getEndPos();
  var map_type: std.posix.system.MAP = .{ .TYPE = .PRIVATE };
  map_type.POPULATE = true;
  const buffer = try std.posix.mmap(
    null,
    size,
    std.posix.PROT.READ,
    map_type,
    file.handle,
    0,
  );
  errdefer std.posix.munmap(buffer);

  return buffer;
}


const Ptype = enum {
  P1, P2, P3, P4, P5, P6,
};

pub fn writePgm(comptime P: Ptype, width: usize, height: usize, pixels: []const u8, filepath: []const u8) !void {
  if (width * height != pixels.len) {
    return error.IncorrectSize;
  }

  // Create file
  const file = try std.fs.cwd().createFile(
    filepath,
    .{ .read = true },
  );
  defer file.close();

  switch (P) {
    .P1 => {
      // Prepare PGM header (https://en.wikipedia.org/wiki/Netpbm)
      var buffer: [255]u8 = [_]u8{ 0 } ** 255;
      const pgmHeader = try std.fmt.bufPrint(&buffer, "P1\n{} {}\n", .{ width, height });
      // Write to file
      try file.writeAll(pgmHeader);
      for (0..height) |j| {
        for (0..width) |i| {
          if (pixels[j * width + i] != 0) {
            try file.writeAll("1 ");
          } else {
            try file.writeAll("0 ");
          }
        }
        try file.writeAll("\n");
      }
    },
    .P6 => {
      // Prepare PGM header (https://en.wikipedia.org/wiki/Netpbm)
      var buffer: [255]u8 = [_]u8{ 0 } ** 255;
      const pgmHeader = try std.fmt.bufPrint(&buffer, "P6\n{} {}\n255\n", .{ width, height });
      // Write to file
      try file.writeAll(pgmHeader);
      const stdout = std.io.getStdOut().writer();
      for (pixels) |p| { try stdout.print("{} ", .{ p }); }
      try file.writeAll(pixels);
    },
    else => return error.NotImplemented,
  }
}

pub fn writePam32(width: usize, height: usize, pixels: []const u8, filepath: []const u8) !void {
  // Create file
  const file = try std.fs.cwd().createFile(
    filepath,
    .{ .read = true },
  );
  defer file.close();
  // Prepare PGM header (https://en.wikipedia.org/wiki/Netpbm)
  var buffer: [255]u8 = [_]u8{ 0 } ** 255;
  const pgmHeader = try std.fmt.bufPrint(&buffer, "P7\nWIDTH {}\nHEIGHT {}\nDEPTH 4\nMAXVAL 255\nTUPLTYPE RGB_ALPHA\nENDHDR\n", .{ width, height });
  // Write to file
  try file.writeAll(pgmHeader);
  try file.writeAll(pixels);
}
