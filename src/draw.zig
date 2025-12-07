const std = @import("std");

const fontfile = @embedFile("dos_8x8_font_white.pbm");

fn asf32(value: anytype) f32 {
  return @as(f32, @floatFromInt(value));
}

const PathCommandType = enum {
  move_to,
  line_to,
  close_path,
};

const PathCommand = union(PathCommandType) {
  move_to: struct { i16, i16 },
  line_to: struct { i16, i16 },
  close_path: void,
};

pub const PixelFormat = enum {
  RGB8,
  RGBA8,
};

// FIXME: Deal with image data other than RGBA
pub const ImageData = struct {
  width: u32,
  height: u32,
  pixel_format: PixelFormat,
  data: []const u8,
};

pub const DrawContext = struct {
  const Self = @This();

  const _a = 0; const _b = 1; const _c = 2; const _d = 3; const _e = 4; const _f = 5;

  width: u32,
  height: u32,
  buffer: []u32,
  strokeStyle: u32 = 0xFFFFFFFF,
  fillStyle: u32 = 0xFFFFFFFF,
  thickness: u32 = 0,
  alpha: bool = false,

  transformMatrix: [6]f32 = [_]f32{ 1, 0, 0, 1, 0, 0 },
  stack: [6]f32 = [_]f32{ 1, 0, 0, 1, 0, 0 },

  path_command_stack: std.ArrayList(PathCommand),

  allocator: std.mem.Allocator,

  pub fn init(allocator: std.mem.Allocator, width: u32, height: u32) !Self {
    const buffer = try allocator.alloc(u32, width * height);
    @memset(buffer, 0);

    return Self{
      .width = width,
      .height = height,
      .buffer = buffer,
      .path_command_stack = std.ArrayList(PathCommand).init(allocator),
      .allocator = allocator,
    };
  }

  pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
    allocator.free(self.buffer);
    self.path_command_stack.deinit();
  }

  // resets (overrides) the current transformation to the identity matrix, and
  // then invokes a transformation described by the arguments of this method.
  // This lets you scale, rotate, translate (move), and skew the context.
  //                       a c e
  // The transformation matrix is described by: [ b d f ]
  //                       0 0 1
  pub fn setTransform(self: *Self, a: f32, b: f32, c: f32, d: f32, e: f32, f: f32) void {
    self.transformMatrix = .{ a, b, c, d, e, f };
  }
  // multiplies the current transformation with the matrix described by the
  // arguments of this method. This lets you scale, rotate, translate (move),
  // and skew the context.
  pub fn transform(self: *Self, a: f32, b: f32, c: f32, d: f32, e: f32, f: f32) void {
    self.transformMatrix[_a] = self.transformMatrix[_a] * a + self.transformMatrix[_c] * b;
    self.transformMatrix[_b] = self.transformMatrix[_b] * a + self.transformMatrix[_d] * b;
    self.transformMatrix[_c] = self.transformMatrix[_a] * c + self.transformMatrix[_c] * d;
    self.transformMatrix[_d] = self.transformMatrix[_b] * c + self.transformMatrix[_d] * d;
    self.transformMatrix[_e] = self.transformMatrix[_a] * e + self.transformMatrix[_c] * f + self.transformMatrix[_e] * 1;
    self.transformMatrix[_f] = self.transformMatrix[_b] * e + self.transformMatrix[_d] * f + self.transformMatrix[_f] * 1;
  }
  // retrieves the current transformation matrix being applied to the context.
  pub fn getTransform(self: Self) [6]f32 {
    return self.transformMatrix;
  }
  // Saves the entire state of the canvas by pushing the current state onto a
  // stack.
  // ⚠️ Only one level of stack for now.
  pub fn save(self: Self) !void {
    self.stack = self.transformMatrix;
  }
  // Restores the most recently saved canvas state by popping the top entry
  // in the drawing state stack. If there is no saved state, this method does
  // nothing.
  pub fn restore(self: Self) void {
    self.transformMatrix = self.stack;
  }
  // resets the rendering context to its default state, allowing it to be
  // reused for drawing something else without having to explicitly reset all
  // the properties.
  pub fn reset(self: Self) void {
    self.transformMatrix = .{ 1, 0, 0, 1, 0, 0 };
  }
  // adds a translation transformation to the current matrix.
  pub fn translate(self: Self, x: f32, y: f32) void {
    self.transformMatrix[_e] += x;
    self.transformMatrix[_f] += y;
  }
  // erases the pixels in a rectangular area by setting them to transparent
  // black.
  // ⚠️ Operates in buffer space (do not take the transformation matrix into account)
  pub fn clearRect(self: Self, x: i16, y: i16, width: i16, height: i16) void {
    if (x == 0 and y == 0 and width == self.width and height == self.height) {
      @memset(&self.buffer, 0xFF000000);
    } else {
      @panic("clearReact on sizes different from the canvas is not yet implemented");
    }
  }
  // fill a rectangle with color
  pub fn fillRect(self: Self, x: i16, y: i16, width: i16, height: i16) void {
    const xc, const yc = self.transformPoint(asf32(x), asf32(y));
    const widthc, const heightc = self.transformVector(asf32(width), asf32(height));
    self.fillPolygon(&.{
      [2]f32{ xc, yc },
      [2]f32{ xc + widthc, yc },
      [2]f32{ xc + widthc, yc + heightc },
      [2]f32{ xc, yc + heightc },
    }, self.fillStyle);
  }
  // Renders a rectangle with a starting point is at (x, y) and whose size is
  // specified by width and height.
  pub fn strokeRect(self: Self, x: i16, y: i16, width: i16, height: i16) void {
    const xc, const yc = self.transformPoint(asf32(x), asf32(y));
    const widthc, const heightc = self.transformVector(asf32(width), asf32(height));

    const xi: i16 = @intFromFloat(@round(xc));
    const yi: i16 = @intFromFloat(@round(yc));
    const widthi: i16 = @intFromFloat(@round(widthc));
    const heighti: i16 = @intFromFloat(@round(heightc));

    self.line(xi, yi, xi + widthi, yi);
    self.line(xi + widthi, yi, xi + widthi, yi + heighti);
    self.line(xi + widthi, yi + heighti, xi, yi + heighti);
    self.line(xi, yi + heighti, xi, yi);
  }
  // Renders a filled circle with a center at (x, y) and a radius
  pub fn fillCircle(self: Self, x: i16, y: i16, radius: i16) void {
    // FIXME: This naive implementation does not take into account the skew than
    // can be introduced by the transformation matrix.
    std.debug.assert(self.transformMatrix[_c] == 0);
    std.debug.assert(self.transformMatrix[_b] == 0);
    const xc, const yc = self.transformPoint(asf32(x), asf32(y));
    const radiusc, _ = self.transformVector(asf32(radius), asf32(0));
    const xi: i16 = @intFromFloat(@round(xc));
    const yi: i16 = @intFromFloat(@round(yc));
    const radiusi: i16 = @intFromFloat(@round(radiusc));
    const radiusi_squared: i32 = @as(i32, @intCast(radiusi)) * radiusi;
    var j = yi - radiusi;
    while (j <= yi + radiusi) {
      var i = xi - radiusi;
      while (i <= xi + radiusi) {
        const dx: i32 = i - xi;
        const dy: i32 = j - yi;
        if (dx * dx + dy * dy <= radiusi_squared) {
          self.internalPlot(i, j, self.fillStyle);
        }
        i += 1;
      }
      j += 1;
    }
  }
  // Draws a line.
  pub fn line(self: Self, startx: i16, starty: i16, endx: i16, endy: i16) void {
    // std.debug.assert(startx >= 0 and starty >= 0 and endx >= 0 and endy >= 0);

    // This is commented because, although debug performances are much better,
    // release performances are worst!
    const builtin = @import("builtin");
    if (builtin.mode == .Debug) {
      // This whole block is an optimization for vertical and horizontal line
      const ux =
        self.transformMatrix[_a] * @as(f32, @floatFromInt(startx)) +
        self.transformMatrix[_c] * @as(f32, @floatFromInt(starty)) + self.transformMatrix[_e];
      var uy =
        self.transformMatrix[_b] * @as(f32, @floatFromInt(startx)) +
        self.transformMatrix[_d] * @as(f32, @floatFromInt(starty)) + self.transformMatrix[_f];
      const vx =
        self.transformMatrix[_a] * @as(f32, @floatFromInt(endx)) +
        self.transformMatrix[_c] * @as(f32, @floatFromInt(endy)) + self.transformMatrix[_e];
      var vy =
        self.transformMatrix[_b] * @as(f32, @floatFromInt(endx)) +
        self.transformMatrix[_d] * @as(f32, @floatFromInt(endy)) + self.transformMatrix[_f];

      if (ux > 0 and uy > 0 and vx > 0 and vy > 0 and
        ux < asf32(self.width) and uy < asf32(self.height) and vx < asf32(self.width) and vy < asf32(self.height)) {
        // If a line is entirely in the canvas
        if (ux == vx) {
          // vertical line
          const x: u16 = @intFromFloat(ux);
          if (uy > vy) std.mem.swap(@TypeOf(uy), &uy, &vy);
          var y: u16 = @intFromFloat(uy);
          while (y < @as(u16, @intFromFloat(vy))) : (y += 1) {
            self.buffer[y * self.width + x] = self.strokeStyle;
          }
          return;
        } else if (uy == vy) {
          // horizontal line
          var startBuffer = @as(u16, @intFromFloat(uy)) * self.width + @as(u16, @intFromFloat(ux));
          var endBuffer = @as(u16, @intFromFloat(vy)) * self.width + @as(u16, @intFromFloat(vx));
          if (startBuffer > endBuffer) std.mem.swap(@TypeOf(startBuffer), &startBuffer, &endBuffer);
          @memset(self.buffer[startBuffer..endBuffer], self.strokeStyle);
          return;
        }
      }
    }
    // Otherwise, we use a general but slow algorithm
    self.drawThickLine(startx, starty, endx, endy);
    // self.drawLineOverlap(startx, starty, endx, endy, 0);
    // self.drawLineWu(startx, starty, endx, endy, 0);
  }
  // Draws a point.
  pub inline fn plot(self: Self, x: i16, y: i16, acolor: u32) void {
    // std.log.debug("plot x {} y {} width {} height {} index {} buffer.len {}", .{
    //  x, y, self.width, self.height,
    //  @as(u16, @bitCast(y)) * width + @as(u16, @bitCast(x)),
    //  buffer.len,
    // });
    const vx = self.transformMatrix[_a] * @as(f32, @floatFromInt(x)) + self.transformMatrix[_c] * @as(f32, @floatFromInt(y)) + self.transformMatrix[_e];
    const vy = self.transformMatrix[_b] * @as(f32, @floatFromInt(x)) + self.transformMatrix[_d] * @as(f32, @floatFromInt(y)) + self.transformMatrix[_f];
    if (vx >= 0 and vx < asf32(self.width) and vy >= 0 and vy < asf32(self.height)) {
      self.buffer[@as(u16, @intFromFloat(vy)) * self.width + @as(u16, @intFromFloat(vx))] = acolor;
    }
  }
  // Writes text at the specified position. x and y specifies the top left
  // corner of the text box to be printed.
  pub fn fillText(self: Self, text: []const u8, x: i16, y: i16) void {
    // @setEvalBranchQuota(10000);
    const fontparams = comptime lbl: {
      // Check we deal with a P1 netpbm file (ASCII text)
      std.debug.assert(fontfile[0] == 'P' and fontfile[1] == '1');
      // Retrieve width and height
      var i = 2;
      while (std.ascii.isWhitespace(fontfile[i])) i += 1;
      var j = i;
      while (!std.ascii.isWhitespace(fontfile[j])) j += 1;
      const fontwidth = try std.fmt.parseInt(usize, fontfile[i..j], 10);
      i = j;
      while (std.ascii.isWhitespace(fontfile[i])) i += 1;
      j = i;
      while (!std.ascii.isWhitespace(fontfile[j])) j += 1;
      const fontheight = try std.fmt.parseInt(usize, fontfile[i..j], 10);
      // Get position of first value
      while (std.ascii.isWhitespace(fontfile[j])) j += 1;
      break :lbl .{ fontwidth, fontheight, j };
    };

    const fontwidth = fontparams[0];
    const fontindex = fontparams[2];

    for (text, 0..) |c, cindex| {
      const cusize = @as(usize, @intCast(c));
      for (0..8) |j| {
        for (0..8) |i| {
          // fontwidth + 1 because of the \n
          if (fontfile[fontindex + j * (fontwidth + 1) + cusize * 8 + i] != '0') {
            self.plot(x + @as(i16, @intCast(i)) + @as(i16, @intCast(cindex * 8)), y + @as(i16, @intCast(j)), self.fillStyle);
          }
        }
      }
    }
  }

  //
  // path functions
  //
  pub fn beginPath(self: Self) void { _ = self; }
  pub fn moveTo(self: *Self, x: i16, y: i16) void {
    self.path_command_stack.append(PathCommand{ .move_to = .{ x, y } }) catch {};
  }
  pub fn lineTo(self: *Self, x: i16, y: i16) void {
    self.path_command_stack.append(PathCommand{ .line_to = .{ x, y } }) catch {};
  }
  pub fn closePath(self: *Self) void {
    self.path_command_stack.append(PathCommand{ .close_path = {} }) catch {};
  }
  pub fn stroke(self: *Self) void {
    var pen: [2]i16 = .{ 0, 0 };
    var first_point: [2]i16 = .{ 0, 0 };
    for (self.path_command_stack.items  ) |command| {
      switch (command) {
        PathCommand.move_to => |position| {
          pen[0] = position[0];
          pen[1] = position[1];
          first_point[0] = pen[0];
          first_point[1] = pen[1];
        },
        PathCommand.line_to => |position| {
          self.line(pen[0], pen[1], position[0], position[1]);
          pen[0] = position[0];
          pen[1] = position[1];
        },
        PathCommand.close_path => {
          self.line(pen[0], pen[1], first_point[0], first_point[1]);
        },
      }
    }
    self.path_command_stack.clearRetainingCapacity();
  }
  pub fn fill(self: *Self) void {
    var pen: [2]i16 = .{ 0, 0 };
    var first_point: [2]i16 = .{ 0, 0 };
    var vertices = std.ArrayList([2]f32).init(self.allocator);
    for (self.path_command_stack.items  ) |command| {
      switch (command) {
        PathCommand.move_to => |position| {
          pen[0] = position[0];
          pen[1] = position[1];
          first_point[0] = pen[0];
          first_point[1] = pen[1];
          vertices.append(self.transformPoint(
            @floatFromInt(pen[0]),
            @floatFromInt(pen[1]))) catch {};
        },
        PathCommand.line_to => |position| {
          pen[0] = position[0];
          pen[1] = position[1];
          vertices.append(self.transformPoint(
            @floatFromInt(pen[0]),
            @floatFromInt(pen[1]))) catch {};
        },
        PathCommand.close_path => {
          pen[0] = first_point[0];
          pen[1] = first_point[1];
          vertices.append(self.transformPoint(
            @floatFromInt(pen[0]),
            @floatFromInt(pen[1]))) catch {};
          self.fillPolygon(vertices.items, self.fillStyle);
          vertices.clearRetainingCapacity();
        },
      }
    }
    if (vertices.items.len > 0) {
      self.fillPolygon(vertices.items, self.fillStyle);
    }
    vertices.deinit();
    self.path_command_stack.clearRetainingCapacity();
  }

  // pub fn putImageData(self: *Self, imageData: ImageData, dx: isize, dy: isize) void {
  //   if (dx <= -@as(i64, @intCast(imageData.width)) or dx >= self.width or
  //       dy <= -@as(i64, @intCast(imageData.height)) or dy >= self.height) {
  //     return;
  //   }

  //   const bounded_dx = @min(self.width - 1, @max(0, dx));
  //   const bounded_dy = @min(self.height - 1, @max(0, dy));

  //   const sx: usize = if (dx >= 0) 0 else @abs(dx);
  //   const sy: usize = if (dy >= 0) 0 else @abs(dy);

  //   const width = if (dx >= 0)
  //     if (bounded_dx + imageData.width < self.width)
  //       imageData.width
  //     else
  //       imageData.width -| ((bounded_dx + imageData.width) % self.width)
  //   else
  //     imageData.width - @abs(dx);
  //   const height = if (dy >= 0)
  //     if (bounded_dy + imageData.height < self.height)
  //       imageData.height
  //     else
  //       imageData.height -| ((bounded_dy + imageData.height) % self.height)
  //   else
  //     imageData.height - @abs(dy);

  //   std.log.debug("dx {} dy {} sx {} sy {} width {} height {}", .{ dx, dy, sx, sy, width, height });
  //   self.putImageData2(imageData, bounded_dx, bounded_dy, sx, sy, width, height);
  // }

  // A little slower be easier to understand
  pub fn putImageData(self: *Self, imageData: ImageData, dx: isize, dy: isize) void {
      if (dx >= self.width or
          dy >= self.height or
          dx <= -@as(i64, @intCast(imageData.width)) or
          dy <= -@as(i64, @intCast(imageData.height))) {
          return;
      }

      const sx: usize = if (dx < 0) @intCast(-dx) else 0;
      const sy: usize = if (dy < 0) @intCast(-dy) else 0;

      const dest_x: usize = if (dx < 0) 0 else @intCast(dx);
      const dest_y: usize = if (dy < 0) 0 else @intCast(dy);

      const width = @min(imageData.width - sx, self.width - dest_x);
      const height = @min(imageData.height - sy, self.height - dest_y);

      self.putImageData2(imageData, dest_x, dest_y, sx, sy, width, height);
  }

  pub fn putImageData2(self: *Self, imageData: ImageData, dx: usize, dy: usize,
    sx: usize, sy: usize, width: usize, height: usize) void {
    const startBuffer: usize = @intCast(dx + dy * self.width);
    if (!self.alpha) {
      const ptr_u32: [*]const u32 = @ptrCast(@alignCast(imageData.data.ptr));
      const rgbadata = ptr_u32[0 .. imageData.data.len / 4];
      for (0..height) |j| {
        const bufferStride = j * self.width;
        const imageDataStride = sx + (sy + j) * imageData.width;
        @memcpy(
          self.buffer[startBuffer + bufferStride..startBuffer + bufferStride + width],
          rgbadata[imageDataStride..imageDataStride + width]
        );
      }
    } else {
      // raw pointers to avoid slice overhead in the loop
      const dest_ptr: [*]u32 = self.buffer.ptr;
      const src_ptr: [*]const u32 = @ptrCast(@alignCast(imageData.data.ptr));
      const start_buffer: usize = @intCast(dx + dy * self.width);
      for (0..height) |j| {
        var buffer_index = start_buffer + j * self.width;
        var image_index = sx + (sy + j) * imageData.width;
        for (0..width) |i| {
          _ = i;
          defer {
            buffer_index += 1;
            image_index += 1;
          }

          const fg_packed = src_ptr[image_index];
          const bg_packed = dest_ptr[buffer_index];
          // We expand to u16 immediately to prevent overflow during multiply
          const fg: @Vector(4, u16) = @intCast(@as(@Vector(4, u8), @bitCast(fg_packed)));
          const bg: @Vector(4, u16) = @intCast(@as(@Vector(4, u8), @bitCast(bg_packed)));
          const alpha_val = fg[3];
          // Optimization: If alpha is 0, skip. If 255, simple copy.
          if (alpha_val == 0) continue;
          if (alpha_val == 255) {
              dest_ptr[buffer_index] = fg_packed;
              continue;
          }

          const a: @Vector(4, u16) = @splat(alpha_val);
          const max: @Vector(4, u16) = @splat(255);
          const inv_a = max - a;
          // (bg * inv_a + fg * a) / 255
          const tmp = bg * inv_a + fg * a;
          const result_16 = (tmp + @as(@Vector(4, u16), @splat(1)) + (tmp >> @as(@Vector(4, u8), @splat(8)))) >> @as(@Vector(4, u8), @splat(8));

          const result_8: @Vector(4, u8) = @intCast(result_16);
          // 7. Store full pixel (Single instruction store)
          dest_ptr[buffer_index] = @bitCast(result_8);
        }
      }
    }
  }

  //
  // Private functions
  //

  /// Draws a point without transformation. To be used by private functions.
  inline fn internalPlot(self: Self, x: i16, y: i16, acolor: u32) void {
    if (x >= 0 and x < self.width and y >= 0 and y < self.height) {
      self.buffer[@as(u32, @intCast(y)) * self.width + @as(u32, @intCast(x))] = acolor;
    }
  }

  /// Fills a polygon defined by `vertices` using the scan-line algorithm.
  fn fillPolygon(self: Self, vertices: []const [2]f32, afillStyle: u32) void {
    const Fns = struct {
      /// Represents a single edge of the polygon, prepared for the scan-line algorithm.
      /// This version stores 'slope' (dy/dx) to match the original Python code.
      const Edge = struct {
        min_y: i16,
        max_y: i16,
        x_at_min_y: f32,
        invslope: f32,
      };

      /// Define a function to calculate the slope (dy/dx) of an edge.
      fn calculateInvSlope(x0: f32, y0: f32, x1: f32, y1: f32) f32 {
        const dx = x0 - x1;
        const dy = y0 - y1;
        if (dy == 0) {
          return std.math.inf(f32); // Horizontal line
        }
        return dx / dy;
      }

      /// Define a function to initialize all edges of the polygon.
      fn initializeEdges(vrtces: []const [2]f32, edges: *std.ArrayListUnmanaged(Edge)) void {
        for (vrtces, 0..) |p1, i| {
          const p2 = vrtces[(i + 1) % vrtces.len];

          // Round y-coordinates to determine integer scan-lines.
          const y0_rounded = @round(p1[1]);
          const y1_rounded = @round(p2[1]);

          const min_y: i16 = @intFromFloat(@min(y0_rounded, y1_rounded));
          const max_y: i16 = @intFromFloat(@max(y0_rounded, y1_rounded));
          const x_at_min_y = if (p1[1] <= p2[1]) p1[0] else p2[0];
          const invslope = calculateInvSlope(p1[0], p1[1], p2[0], p2[1]);

          edges.appendAssumeCapacity(Edge{
            .min_y = min_y,
            .max_y = max_y,
            .x_at_min_y = x_at_min_y,
            .invslope = invslope,
          });
        }
      }

      /// Define a function to initialize the global edge table.
      fn initializeGlobalEdgeTable(edges: []const Edge, global_edge_table: *std.ArrayListUnmanaged(Edge)) void {
        global_edge_table.appendSliceAssumeCapacity(edges);

        // Sort by min_y, then x_at_min_y
        const sort_context = struct {
          fn lessThan(_: void, a: Edge, b: Edge) bool {
            if (a.min_y != b.min_y) {
              return a.min_y < b.min_y;
            }
            return a.x_at_min_y < b.x_at_min_y;
          }
        };
        std.mem.sort(Edge, global_edge_table.items, {}, sort_context.lessThan);

        // We eliminate flat edges (slope == 0)
        var i: usize = 0;
        var len = global_edge_table.items.len;
        while (i < len) {
          if (global_edge_table.items[i].invslope == std.math.inf(f32)) {
            _ = global_edge_table.orderedRemove(i);
            len = global_edge_table.items.len;
          } else {
            i += 1;
          }
        }
      }

      /// Select the edges that will be used for the scanline `scan_line`.
      /// This function uses the original Python's flawed sorting logic.
      fn getActiveEdgeTable(
        source_table: []const Edge,
        scan_line: i32,
        active_edge_table: *std.ArrayListUnmanaged(Edge),
      ) void {
        active_edge_table.clearRetainingCapacity();
        for (source_table) |edge| {
          if (edge.min_y <= scan_line and edge.max_y > scan_line) {
            active_edge_table.appendAssumeCapacity(edge);
          }
        }

        // Sort by the static key (x_at_min_y, -slope), matching the original Python code.
        const sort_context = struct {
          fn lessThan(_: void, a: Edge, b: Edge) bool {
            if (a.x_at_min_y != b.x_at_min_y) {
              return a.x_at_min_y < b.x_at_min_y;
            }
            return a.invslope < b.invslope;
          }
        };
        std.mem.sort(Edge, active_edge_table.items, {}, sort_context.lessThan);
      }
    };

    if (vertices.len < 3) return;

    var edge_buffer: [64]Fns.Edge = undefined; // Preallocate a buffer for edges.
    var edge_list = std.ArrayListUnmanaged(Fns.Edge).initBuffer(&edge_buffer);
    Fns.initializeEdges(vertices, &edge_list);

    var global_edge_buffer: [64]Fns.Edge = undefined; // Preallocate a buffer for edges.
    var global_edge_table = std.ArrayListUnmanaged(Fns.Edge).initBuffer(&global_edge_buffer);
    Fns.initializeGlobalEdgeTable(edge_list.items, &global_edge_table);

    if (global_edge_table.items.len == 0) return;

    // Initialize scan-line and active edge table.
    var scan_line: i16 = global_edge_table.items[0].min_y;
    var active_edge_table = std.ArrayListUnmanaged(Fns.Edge).initBuffer(&edge_buffer);
    Fns.getActiveEdgeTable(global_edge_table.items, scan_line, &active_edge_table);

    // Iterate over each scan-line until active edge table is empty.
    while (active_edge_table.items.len > 0) {
      // Draw pixels between x-values of odd and even parity edge pairs.
      var i: usize = 0;
      while (i + 1 < active_edge_table.items.len) : (i += 2) {
        const even_edge = active_edge_table.items[i];
        const odd_edge = active_edge_table.items[i + 1];

        const x_start_f = asf32(scan_line - even_edge.min_y) * even_edge.invslope + even_edge.x_at_min_y;
        const x_end_f = asf32(scan_line - odd_edge.min_y) * odd_edge.invslope + odd_edge.x_at_min_y;

        var x_start = @round(x_start_f);
        const x_end = @round(x_end_f);

        while (x_start < x_end) {
          self.internalPlot(@intFromFloat(x_start), scan_line, afillStyle);
          x_start += 1;
        }
      }

      // Update scan-line and rebuild active edge table for the next line.
      scan_line += 1;
      Fns.getActiveEdgeTable(global_edge_table.items, scan_line, &active_edge_table);
    }
  }

  /// Transform a point according to the currently set transform matrix in the context.
  fn transformPoint(self: Self, x: f32, y: f32) struct { f32, f32 } {
    return .{
      self.transformMatrix[_a] * x + self.transformMatrix[_c] * y + self.transformMatrix[_e],
      self.transformMatrix[_b] * x + self.transformMatrix[_d] * y + self.transformMatrix[_f],
    };
  }

  /// Transform a vector according to the currently set transform matrix in the context.
  fn transformVector(self: Self, x: f32, y: f32) struct { f32, f32 } {
    return .{
      self.transformMatrix[_a] * x + self.transformMatrix[_c] * y,
      self.transformMatrix[_b] * x + self.transformMatrix[_d] * y,
    };
  }

  // Modified Bresenham draw(line) with optional overlap. Required for drawThickLine().
  // Overlap draws additional pixel when changing minor direction. For standard bresenham overlap, choose LINE_OVERLAP_NONE (0).
  //
  // Sample line:
  //
  //  00+
  //   -0000+
  //     -0000+
  //       -00
  //
  // 0 pixels are drawn for normal line without any overlap LINE_OVERLAP_NONE
  // + pixels are drawn if LINE_OVERLAP_MAJOR
  // - pixels are drawn if LINE_OVERLAP_MINOR

  const Overlap = enum(u8) {
    LINE_OVERLAP_NONE = 0,
    LINE_OVERLAP_MINOR = 1,
    LINE_OVERLAP_MAJOR = 2,
    LINE_OVERLAP_BOTH = 3,
  };

  fn drawLineOverlap(self: Self, pstartx: i16, pstarty: i16, endx: i16, endy: i16, aOverlap: u8) void {
    var tStepX: i16 = 0;
    var tStepY: i16 = 0;
    var tDeltaXTimes2: i16 = 0;
    var tDeltaYTimes2: i16 = 0;
    var tError: i16 = 0;
    var startx = pstartx;
    var starty = pstarty;
    // calculate direction
    var tDeltaX = endx - startx;
    var tDeltaY = endy - starty;
    if (tDeltaX < 0) {
      tDeltaX = -tDeltaX;
      tStepX = -1;
    } else {
      tStepX = 1;
    }
    if (tDeltaY < 0) {
      tDeltaY = -tDeltaY;
      tStepY = -1;
    } else {
      tStepY = 1;
    }
    tDeltaXTimes2 = tDeltaX << 1;
    tDeltaYTimes2 = tDeltaY << 1;
    // draw start pixel
    self.plot(startx, starty, self.strokeStyle);
    if (tDeltaX > tDeltaY) {
      // start value represents a half step in Y direction
      tError = tDeltaYTimes2 - tDeltaX;
      while (startx != endx) {
        // step in main direction
        startx += tStepX;
        if (tError >= 0) {
          if (aOverlap & @intFromEnum(Overlap.LINE_OVERLAP_MAJOR) != 0) {
            // draw pixel in main direction before changing
            self.plot(startx, starty, self.strokeStyle);
          }
          // change Y
          starty += tStepY;
          if (aOverlap & @intFromEnum(Overlap.LINE_OVERLAP_MINOR) != 0) {
            // draw pixel in minor direction before changing
            self.plot(startx - tStepX, starty, self.strokeStyle);
          }
          tError -= tDeltaXTimes2;
        }
        tError += tDeltaYTimes2;
        self.plot(startx, starty, self.strokeStyle);
      }
    } else {
      tError = tDeltaXTimes2 - tDeltaY;
      while (starty != endy) {
        starty += tStepY;
        if (tError >= 0) {
          if (aOverlap & @intFromEnum(Overlap.LINE_OVERLAP_MAJOR) != 0) {
            // draw pixel in main direction before changing
            self.plot(startx, starty, self.strokeStyle);
          }
          startx += tStepX;
          if (aOverlap & @intFromEnum(Overlap.LINE_OVERLAP_MINOR) != 0) {
            // draw pixel in minor direction before changing
            self.plot(startx, starty - tStepY, self.strokeStyle);
          }
          tError -= tDeltaYTimes2;
        }
        tError += tDeltaXTimes2;
        self.plot(startx, starty, self.strokeStyle);
      }
    }
  }

  // fractional part of x
  fn fpart(x: f32) f32 {
    return x - @floor(x);
  }

  fn rfpart(x: f32) f32 {
    return 1 - fpart(x);
  }

  pub fn shadeColor(R: u8, G: u8, B: u8) u32 {
    return R | (@as(u32, G) << 8) | (@as(u32, B) << 16) | 0xFF000000; // Always full alpha
  }

  fn drawLineWu(self: Self, px0: i16, py0: i16, px1: i16, py1: i16, unused: u8) void {
    std.log.debug("drawLineWu {} {} {} {}", .{ px0, py0, px1, py1 });
    _ = unused;
    var x0: f32 = @floatFromInt(px0);
    var y0: f32 = @floatFromInt(py0);
    var x1: f32 = @floatFromInt(px1);
    var y1: f32 = @floatFromInt(py1);
    const R: f32 = @floatFromInt(self.strokeStyle & 0x000000FF);
    const G: f32 = @floatFromInt((self.strokeStyle & 0x0000FF00) >> 8);
    const B: f32 = @floatFromInt((self.strokeStyle & 0x00FF0000) >> 16);

    const steep = @abs(y1 - y0) > @abs(x1 - x0);

    if (steep) {
      std.mem.swap(f32, &x0, &y0);
      std.mem.swap(f32, &x1, &y1);
      std.log.debug("drawLineWu2 {d} {d} {d} {d}", .{ x0, y0, x1, y1 });
    }
    if (x0 > x1) {
      std.mem.swap(f32, &x0, &x1);
      std.mem.swap(f32, &y0, &y1);
      std.log.debug("drawLineWu3 {d} {d} {d} {d}", .{ x0, y0, x1, y1 });
    }

    const dx = x1 - x0;
    const dy = y1 - y0;

    var gradient: f32 = 1.0;
    if (dx != 0.0) {
      gradient = dy / dx;
    }

    // handle first endpoint
    var xend = @round(x0);
    var yend = y0 + gradient * (xend - x0);
    var xgap = rfpart(x0 + 0.5);
    const xpxl1: i16 = @intFromFloat(xend); // this will be used in the main loop
    const ypxl1: i16 = @intFromFloat(@floor(yend));
    if (steep) {
      self.plot(ypxl1  , xpxl1  , shadeColor(R, G, B, rfpart(yend) * xgap));
      self.plot(ypxl1 + 1, xpxl1  , shadeColor(R, G, B, fpart(yend) * xgap));
    } else {
      self.plot(xpxl1  , ypxl1  , shadeColor(R, G, B, rfpart(yend) * xgap));
      self.plot(xpxl1  , ypxl1 + 1, shadeColor(R, G, B, fpart(yend) * xgap));
    }
    var intery = yend + gradient; // first y-intersection for the main loop

    // handle second endpoint
    std.log.debug("x1 {d} y1 {d}", .{ x1, y1 });
    xend = @round(x1);
    yend = y1 + gradient * (xend - x1);
    std.log.debug("xend {d} yend {d}", .{ xend, yend });
    xgap = fpart(x1 + 0.5);
    const xpxl2: i16 = @intFromFloat(xend); //this will be used in the main loop
    const ypxl2: i16 = @intFromFloat(@floor(yend));
    std.log.debug("xpxl2 {} ypxl2 {}", .{ xpxl2, ypxl2 });
    if (steep) {
      self.plot(ypxl2  , xpxl2  , shadeColor(R, G, B, rfpart(yend) * xgap));
      self.plot(ypxl2 + 1, xpxl2  , shadeColor(R, G, B, fpart(yend) * xgap));
    } else {
      self.plot(xpxl2  , ypxl2  , shadeColor(R, G, B, rfpart(yend) * xgap));
      self.plot(xpxl2  , ypxl2 + 1, shadeColor(R, G, B, fpart(yend) * xgap));
    }

    // main loop
    if (steep) {
      var i = xpxl1 + 1;
      while (i < xpxl2 - 1) {
        self.plot(@as(i16, @intFromFloat(@floor(intery)))  , i, shadeColor(R, G, B, rfpart(intery)));
        self.plot(@as(i16, @intFromFloat(@floor(intery))) + 1, i, shadeColor(R, G, B, fpart(intery)));
        intery = intery + gradient;
        i += 1;
      }
    } else {
      var i = xpxl1 + 1;
      while (i < xpxl2 - 1) {
        self.plot(i, @as(i16, @intFromFloat(@floor(intery)))  , shadeColor(R, G, B, rfpart(intery)));
        self.plot(i, @as(i16, @intFromFloat(@floor(intery))) + 1, shadeColor(R, G, B, fpart(intery)));
        intery = intery + gradient;
        i += 1;
      }
    }
  }

  //
  // The same as before, but no clipping to display range, some pixel are drawn twice (because of using LINE_OVERLAP_BOTH)
  // and direction of thickness changes for each octant (except for LINE_THICKNESS_MIDDLE and thickness value is odd)
  // thicknessMode can be LINE_THICKNESS_MIDDLE or any other value
  //
  fn drawThickLine(self: Self, pstartx: i16, pstarty: i16, pendx: i16, pendy: i16) void {
    var tStepX: i16 = 0;
    var tStepY: i16 = 0;
    var tDeltaXTimes2: i16 = 0;
    var tDeltaYTimes2: i16 = 0;
    var tError: i16 = 0;
    var startx = pstartx;
    var starty = pstarty;
    var endx = pendx;
    var endy = pendy;

    var tDeltaY = startx - endx;
    var tDeltaX = endy - startx;
    // mirror 4 quadrants to one and adjust deltas and stepping direction
    if (tDeltaX < 0) {
      tDeltaX = -tDeltaX;
      tStepX = -1;
    } else {
      tStepX = 1;
    }
    if (tDeltaY < 0) {
      tDeltaY = -tDeltaY;
      tStepY = -1;
    } else {
      tStepY = 1;
    }
    tDeltaXTimes2 = tDeltaX << 1;
    tDeltaYTimes2 = tDeltaY << 1;
    var tOverlap: Overlap = Overlap.LINE_OVERLAP_NONE;
    // which octant are we now
    if (tDeltaX > tDeltaY) {
      // if (we want to draw the original coordinate in the middle of the thick line)
      {
        // adjust draw start point
        tError = tDeltaYTimes2 - tDeltaX;
        var i = self.thickness / 2;
        while (i > 0) {
          // change X (main direction here)
          startx -= tStepX;
          endx -= tStepX;
          if (tError >= 0) {
            // change Y
            starty -= tStepY;
            endy -= tStepY;
            tError -= tDeltaXTimes2;
          }
          tError += tDeltaYTimes2;
          i -= 1;
        }
      }
      self.drawLineOverlap(startx, starty, endx, endy, @intFromEnum(tOverlap));
      // drawLineWu(startx, starty, endx, endy, @intFromEnum(tOverlap));
      // draw thickness lines
      tError = tDeltaYTimes2 - tDeltaX;
      var i = self.thickness;
      while (i > 1) {
        // change X (main direction here)
        startx += tStepX;
        endx += tStepX;
        tOverlap = Overlap.LINE_OVERLAP_NONE;
        if (tError >= 0) {
          // change Y
          startx += tStepY;
          endy += tStepY;
          tError -= tDeltaXTimes2;
          tOverlap = Overlap.LINE_OVERLAP_BOTH;
        }
        tError += tDeltaYTimes2;
        self.drawLineOverlap(startx, starty, endx, endy, @intFromEnum(tOverlap));
        // drawLineWu(startx, starty, endx, endy, @intFromEnum(tOverlap));
        i -= 1;
      }
    } else {
      // if (we want to draw the original coordinate in the middle of the thick line)
      {
        tError = tDeltaXTimes2 - tDeltaY;
        var i = self.thickness / 2;
        while (i > 0) {
          starty -= tStepY;
          endy -= tStepY;
          if (tError >= 0) {
            startx -= tStepX;
            endx -= tStepX;
            tError -= tDeltaYTimes2;
          }
          tError += tDeltaXTimes2;
          i -= 1;
        }
      }
      self.drawLineOverlap(startx, starty, endx, endy, @intFromEnum(tOverlap));
      // drawLineWu(startx, starty, endx, endy, @intFromEnum(tOverlap));
      tError = tDeltaXTimes2 - tDeltaY;
      var i = self.thickness;
      while (i > 1) {
        starty += tStepY;
        endy += tStepY;
        tOverlap = Overlap.LINE_OVERLAP_NONE;
        if (tError >= 0) {
          startx += tStepX;
          endx += tStepX;
          tError -= tDeltaYTimes2;
          tOverlap = Overlap.LINE_OVERLAP_BOTH;
        }
        tError += tDeltaXTimes2;
        self.drawLineOverlap(startx, starty, endx, endy, @intFromEnum(tOverlap));
        // drawLineWu(startx, starty, endx, endy, @intFromEnum(tOverlap));
        i -= 1;
      }
    }
  }
};

fn printBuffer(ctx: DrawContext) !void {
  var buffer: [1024]u8 = undefined;
  var b: usize = 0;
  b += (try std.fmt.bufPrint(buffer[b..], "\n", .{})).len;
  for (0..ctx.height) |j| {
    for (0..ctx.width) |i| {
      b += (try std.fmt.bufPrint(buffer[b..], "{} ", .{ ctx.buffer[j * ctx.width + i] })).len;
    }
    b += (try std.fmt.bufPrint(buffer[b..], "\n", .{})).len;
  }
  std.log.debug("{s}", .{ buffer[0..b] });
}

test "fillPolygon" {
  std.testing.log_level = .debug;
  const allocator = std.testing.allocator;

  var vertices = [_][2]f32{ [2]f32{ 1.0, 1.0 }, [2]f32{ 3.0, 1.0 }, [2]f32{ 3.0, 3.0 }, [2]f32{ 1.0, 3.0 } };
  var ctx = try DrawContext.init(allocator, 5, 5);
  defer ctx.deinit(allocator);
  ctx.fillPolygon(&vertices, 1);
  try printBuffer(ctx);
  // FIXME: We shouldn't paint the edge of the polygon
  try std.testing.expectEqualSlices(u32, &.{
    // 0, 0, 0, 0, 0,
    // 0, 0, 0, 0, 0,
    // 0, 0, 1, 0, 0,
    // 0, 0, 0, 0, 0,
    // 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 1, 1, 0, 0,
    0, 1, 1, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
  }, ctx.buffer);
}

test "fillCircle" {
  std.testing.log_level = .debug;
  const allocator = std.testing.allocator;

  var ctx = try DrawContext.init(allocator, 5, 5);
  defer ctx.deinit(allocator);
  ctx.fillStyle = 1;
  ctx.fillCircle(2, 2, 2);
  try printBuffer(ctx);
  try std.testing.expectEqualSlices(u32, &.{
    0, 0, 1, 0, 0,
    0, 1, 1, 1, 0,
    1, 1, 1, 1, 1,
    0, 1, 1, 1, 0,
    0, 0, 1, 0, 0,
  }, ctx.buffer);
}
