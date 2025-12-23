// zig test -freference-trace=8 src/draw.zig
const std = @import("std");

const fontfile = @embedFile("dos_8x8_font_white.pbm");

fn asf32(value: anytype) f32 {
  return @as(f32, @floatFromInt(value));
}

const PathCommandType = enum {
  move_to,
  line_to,
  quad_to,
  close_path,
};

const PathCommand = union(PathCommandType) {
  move_to: struct { i16, i16 },
  line_to: struct { i16, i16 },
  quad_to: struct { cpx: i16, cpy: i16, x: i16, y: i16 },
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
  premultipliedAlpha: bool = false,
};

const InterpolationType = enum {
  nearest,
  linear,
  cubic,
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
  line_dash_segments: []const f32 = &[_]f32{},
  alpha: bool = false,
  // For now, always considered false
  imageSmoothingEnabled: bool = false,
  imageSmoothingQuality: []const u8 = "medium",

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
  //                                              a c e
  // The transformation matrix is described by: [ b d f ]
  //                                              0 0 1
  pub fn setTransform(self: *Self, a: f32, b: f32, c: f32, d: f32, e: f32, f: f32) void {
    self.transformMatrix = .{ a, b, c, d, e, f };
  }
  // multiplies the current transformation with the matrix described by the
  // arguments of this method. This lets you scale, rotate, translate (move),
  // and skew the context.
  pub fn transform(self: *Self, a: f32, b: f32, c: f32, d: f32, e: f32, f: f32) void {
    const tmp_a = self.transformMatrix[_a] * a + self.transformMatrix[_c] * b;
    const tmp_b = self.transformMatrix[_b] * a + self.transformMatrix[_d] * b;
    const tmp_c = self.transformMatrix[_a] * c + self.transformMatrix[_c] * d;
    const tmp_d = self.transformMatrix[_b] * c + self.transformMatrix[_d] * d;
    const tmp_e = self.transformMatrix[_a] * e + self.transformMatrix[_c] * f + self.transformMatrix[_e] * 1;
    const tmp_f = self.transformMatrix[_b] * e + self.transformMatrix[_d] * f + self.transformMatrix[_f] * 1;
    self.transformMatrix[_a] = tmp_a;
    self.transformMatrix[_b] = tmp_b;
    self.transformMatrix[_c] = tmp_c;
    self.transformMatrix[_d] = tmp_d;
    self.transformMatrix[_e] = tmp_e;
    self.transformMatrix[_f] = tmp_f;
  }
  // Apply a rotation by transformaing the matrix.
  // angle: The rotation angle, clockwise in radians. You can use
  // degree * Math.PI / 180 to calculate a radian from a degree.
  pub fn rotate(self: *Self, angle: f32) void {
    const cos = std.math.cos(angle);
    const sin = std.math.sin(angle);
    self.transform(cos, sin, -sin, cos, 0, 0);
  }
  // Apply a scaling matrix
  // x: Scaling factor in the horizontal direction. A negative value flips
  // pixels across the vertical axis. A value of 1 results in no horizontal
  // scaling.
  // y: Scaling factor in the vertical direction. A negative value
  // flips pixels across the horizontal axis. A value of 1 results in no
  // vertical scaling.
  pub fn scale(self: *Self, x: f32, y: f32) void {
    self.transform(x, 0, 0, y, 0, 0);
  }
  // retrieves the current transformation matrix being applied to the context.
  pub fn getTransform(self: Self) [6]f32 {
    return self.transformMatrix;
  }
  // Get the inverse of the transform (canvas → user)
  //
  // det = a*d - b*c
  // if det == 0: return   // non-invertible
  //
  // invDet = 1 / det
  //
  // a' =  d * invDet
  // b' = -b * invDet
  // c' = -c * invDet
  // d' =  a * invDet
  // e' = (c*f - d*e) * invDet
  // f' = (b*e - a*f) * invDet
  //
  // Inverse:
  // | a' c' e' |
  // | b' d' f' |
  // | 0  0  1  |
  fn getInverse(self: Self, inverse: *[6]f32) !void {
    const t = self.transformMatrix;
    const det = t[_a] * t[_d] - t[_b] * t[_c];
    if (det < 0.000001) return error.NotInvertible;

    const inv_det = 1 / det;

    inverse[_a] =  t[_d] * inv_det;
    inverse[_b] = -t[_b] * inv_det;
    inverse[_c] = -t[_c] * inv_det;
    inverse[_d] =  t[_a] * inv_det;
    inverse[_e] = (t[_c] * t[_f] - t[_d] * t[_e]) * inv_det;
    inverse[_f] = (t[_b] * t[_e] - t[_a] * t[_f]) * inv_det;
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
  pub fn translate(self: *Self, x: f32, y: f32) void {
    self.transform(1, 0, 0, 1, x, y);
  }
  // erases the pixels in a rectangular area by setting them to transparent
  // black.
  // ⚠️ Operates in buffer space (do not take the transformation matrix into account)
  pub fn clearRect(self: Self, x: i16, y: i16, width: usize, height: usize) void {
    if (x == 0 and y == 0 and width == self.width and height == self.height) {
      @memset(self.buffer, 0xFF000000);
    } else {
      @panic("clearReact on sizes different from the canvas is not yet implemented");
    }
  }
  // fill a rectangle with color
  pub fn fillRect(self: Self, x: i16, y: i16, width: i16, height: i16) void {
    const xc, const yc = self.transformPoint(asf32(x), asf32(y));
    const widthc, const heightc = self.transformVector(asf32(width), asf32(height));
    self.fillPolygon(self.allocator, &.{
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
    const ux: i16 = @intFromFloat(
      self.transformMatrix[_a] * @as(f32, @floatFromInt(startx)) +
      self.transformMatrix[_c] * @as(f32, @floatFromInt(starty)) + self.transformMatrix[_e]);
    const uy: i16 = @intFromFloat(
      self.transformMatrix[_b] * @as(f32, @floatFromInt(startx)) +
      self.transformMatrix[_d] * @as(f32, @floatFromInt(starty)) + self.transformMatrix[_f]);
    const vx: i16 = @intFromFloat(
      self.transformMatrix[_a] * @as(f32, @floatFromInt(endx)) +
      self.transformMatrix[_c] * @as(f32, @floatFromInt(endy)) + self.transformMatrix[_e]);
    const vy: i16 = @intFromFloat(
      self.transformMatrix[_b] * @as(f32, @floatFromInt(endx)) +
      self.transformMatrix[_d] * @as(f32, @floatFromInt(endy)) + self.transformMatrix[_f]);
    self.drawThickLine(ux, uy, vx, vy);
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

  pub fn setLineDash(self: *Self, segments: []const f32) void {
    self.line_dash_segments = segments;
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
  pub fn quadraticCurveTo(self: *Self, cpx: i16, cpy: i16, x: i16, y: i16) void {
    self.path_command_stack.append(PathCommand{ .quad_to = .{ .x = x, .y = y, .cpx= cpx, .cpy = cpy } }) catch {};
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
        PathCommand.quad_to => |parameters| {
          const start_point: @Vector(2, f32) = @floatFromInt(@as(@Vector(2, i16), pen));
          const control_point = @as(@Vector(2, f32), @floatFromInt(@Vector(2, i16){ parameters.cpx, parameters.cpy }));
          const end_point = @as(@Vector(2, f32), @floatFromInt(@Vector(2, i16){ parameters.x, parameters.y }));
          var current_position = start_point;
          const step = 10;
          for (0..step) |t_int| {
            const t: f32 = asf32(t_int) / (step - 1);
            const first_segment = lerp(
              start_point,
              control_point,
              t);
            const second_segment = lerp(
              control_point,
              end_point,
              t);
            const curve = lerp(first_segment, second_segment, t);
            self.line(
              @intFromFloat(@round(current_position[0])),
              @intFromFloat(@round(current_position[1])),
              @intFromFloat(@round(curve[0])),
              @intFromFloat(@round(curve[1])));
            current_position = curve;
          }
          pen[0] = parameters.x;
          pen[1] = parameters.y;
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
    for (self.path_command_stack.items) |command| {
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
        PathCommand.quad_to => |parameters| {
          const start_point: @Vector(2, f32) = @floatFromInt(@as(@Vector(2, i16), pen));
          const control_point = @as(@Vector(2, f32), @floatFromInt(@Vector(2, i16){ parameters.cpx, parameters.cpy }));
          const end_point = @as(@Vector(2, f32), @floatFromInt(@Vector(2, i16){ parameters.x, parameters.y }));
          vertices.append(self.transformPoint(start_point[0], start_point[1])) catch {};
          const step = 10;
          for (0..step) |t_int| {
            const t: f32 = asf32(t_int) / (step - 1);
            const first_segment = lerp(
              start_point,
              control_point,
              t);
            const second_segment = lerp(
              control_point,
              end_point,
              t);
            const curve = lerp(first_segment, second_segment, t);
            vertices.append(self.transformPoint(curve[0], curve[1])) catch {};
          }
          pen[0] = parameters.x;
          pen[1] = parameters.y;
        },
        PathCommand.close_path => {
          pen[0] = first_point[0];
          pen[1] = first_point[1];
          vertices.append(self.transformPoint(
            @floatFromInt(pen[0]),
            @floatFromInt(pen[1]))) catch {};
          self.fillPolygon(self.allocator, vertices.items, self.fillStyle);
          vertices.clearRetainingCapacity();
        },
      }
    }
    if (vertices.items.len > 0) {
      self.fillPolygon(self.allocator, vertices.items, self.fillStyle);
    }
    vertices.deinit();
    self.path_command_stack.clearRetainingCapacity();
  }

  // A little faster but more complicated implementation
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

  // A little slower but easier to understand
  // Made with Gemini Pro 3
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
          // // We expand to u16 immediately to prevent overflow during multiply
          // const fg: @Vector(4, u16) = @intCast(@as(@Vector(4, u8), @bitCast(fg_packed)));
          // const bg: @Vector(4, u16) = @intCast(@as(@Vector(4, u8), @bitCast(bg_packed)));
          // const alpha_val = fg[3];
          // // Optimization: If alpha is 0, skip. If 255, simple copy.
          // if (alpha_val == 0) continue;
          // if (alpha_val == 255) {
          //     dest_ptr[buffer_index] = fg_packed;
          //     continue;
          // }

          // const a: @Vector(4, u16) = @splat(alpha_val);
          // const max: @Vector(4, u16) = @splat(255);
          // const inv_a = max - a;
          // // (bg * inv_a + fg * a) / 255
          // const tmp = bg * inv_a + fg * a;
          // const result_16 = (tmp + @as(@Vector(4, u16), @splat(1)) + (tmp >> @as(@Vector(4, u8), @splat(8)))) >> @as(@Vector(4, u8), @splat(8));

          // const result_8: @Vector(4, u8) = @intCast(result_16);
          const result_8 = blend(bg_packed, fg_packed);
          // 7. Store full pixel (Single instruction store)
          dest_ptr[buffer_index] = result_8;
        }
      }
    }
  }

  // drawImage Rasterization Algorithm (Affine Canvas)
  //
  // Spaces:
  // - Image space: source pixels
  // - User space: drawImage(dx,dy,dw,dh)
  // - Canvas space: device pixels
  //
  // Current transform M (user → canvas):
  // | a c e |
  // | b d f |
  // | 0 0 1 |
  //
  // 1) Invert transform (canvas → user)
  //
  // see getInverse
  //
  // 2) Transform destination rect corners to canvas space
  //
  // P0 = M * (dx, dy)
  // P1 = M * (dx + dw, dy)
  // P2 = M * (dx, dy + dh)
  // P3 = M * (dx + dw, dy + dh)
  //
  // 3) Compute canvas-space bounding box of P0–P3
  //    Clip bounds to canvas dimensions
  //
  // 4) For each canvas pixel (X,Y) in bounding box:
  //
  // cx = X + 0.5
  // cy = Y + 0.5
  //
  // // canvas → user
  // ux = a' * cx + c' * cy + e'
  // uy = b' * cx + d' * cy + f'
  //
  // // reject pixels outside destination rect
  // if ux < dx or ux >= dx + dw or
  //    uy < dy or uy >= dy + dh:
  //     continue
  //
  // // user → image
  // tx = (ux - dx) / dw
  // ty = (uy - dy) / dh
  //
  // srcX = sx + tx * sWidth
  // srcY = sy + ty * sHeight
  //
  // if srcX < sx or srcX >= sx + sWidth or
  //    srcY < sy or srcY >= sy + sHeight:
  //     continue
  //
  // // sample source image
  // ix = floor(srcX)
  // iy = floor(srcY)
  // color = image[ix, iy]
  //
  // // write result
  // canvas[X, Y] = color
  //
  // Notes:
  // - Backward mapping prevents holes
  // - Per-pixel inverse avoids FP accumulation
  // - Works for scale/rotate/skew/translate
  // - One matrix inversion per draw call
  pub fn drawImage(self: *Self, image: ImageData, dx: i32, dy: i32) void {
    self.drawImage3(image, 0, 0, image.width, image.height, dx, dy, image.width, image.height);
  }
  pub fn drawImage2(self: *Self, image: ImageData, dx: i32, dy: i32, dWidth: usize, dHeight: usize) void {
    self.drawImage3(image, 0, 0, image.width, image.height, dx, dy, dWidth, dHeight);
  }
  // We generate all the combination of the compile time parameters to generate
  // all the possible functions.
  pub fn drawImage3(self: *Self, image: ImageData, sx: i32, sy: i32, sWidth: usize, sHeight: usize,
    dx: i32, dy: i32, dWidth: usize, dHeight: usize) void {

    const no_stretching = sWidth == dWidth and sHeight == dHeight;
    const runtime_interp = self.getInterpolationType();

    // 1. Iterate over possible boolean values at compile-time
    inline for (.{ false, true }) |no_stretching_param| {
      // 2. Check if the runtime value matches the current compile-time value
      if (no_stretching == no_stretching_param) {
        inline for (.{ false, true }) |canvas_alpha| {
          if (self.alpha == canvas_alpha) {
            inline for (.{ false, true }) |apply_premultiplied_alpha| {
              // We apply the premultiplication only if it is NOT already apply in the image
              if (image.premultipliedAlpha != apply_premultiplied_alpha) {
                // 3. Use 'inline else' to unroll the enum switch
                switch (runtime_interp) {
                  inline else => |it| {
                    self.innerDrawImage3(
                      no_stretching_param, canvas_alpha, apply_premultiplied_alpha, it,
                      image, sx, sy, sWidth, sHeight, dx, dy, dWidth, dHeight
                    );
                    return;
                  },
                }
              }
            }
          }
        }
      }
    }
  }

  pub fn innerDrawImage3(self: *Self, comptime no_stretching: bool, comptime alpha: bool,
    comptime premultiplyAlpha: bool, comptime interpolation_type: InterpolationType,
    image: ImageData, sx: i32, sy: i32, sWidth: usize, sHeight: usize,
    dx: i32, dy: i32, dWidth: usize, dHeight: usize) void {
    const fsx: f32 = @floatFromInt(sx);
    const fsy: f32 = @floatFromInt(sy);
    const fsWidth: f32 = @floatFromInt(sWidth);
    const fsHeight: f32 = @floatFromInt(sHeight);
    const fdx: f32 = @floatFromInt(dx);
    const fdy: f32 = @floatFromInt(dy);
    const fdWidth: f32 = @floatFromInt(dWidth);
    const fdHeight: f32 = @floatFromInt(dHeight);
    const fcanvas_width: f32 = @floatFromInt(self.width);
    const fcanvas_height: f32 = @floatFromInt(self.height);

    // Inverse the matrix transform
    var inverse: [6]f32 = .{ 0 } ** 6;
    self.getInverse(&inverse) catch @panic("non invertible transform");
    // Compute the destination parallelogram
    const p0x, const p0y = self.transformPoint(fdx, fdy);
    const p1x, const p1y = self.transformPoint(fdx + fdWidth, fdy);
    const p2x, const p2y = self.transformPoint(fdx, fdy + fdHeight);
    const p3x, const p3y = self.transformPoint(fdx + fdWidth, fdy + fdHeight);
    // Compute the bounding box taking into account the canvas dimentions
    const minx: usize = @intFromFloat(@round(@max(0, @min(fcanvas_width, @reduce(.Min, @Vector(4, f32){ p0x, p1x, p2x, p3x })))));
    const maxx: usize = @intFromFloat(@round(@max(0, @min(fcanvas_width, @reduce(.Max, @Vector(4, f32){ p0x, p1x, p2x, p3x })))));
    const miny: usize = @intFromFloat(@round(@max(0, @min(fcanvas_height, @reduce(.Min, @Vector(4, f32){ p0y, p1y, p2y, p3y })))));
    const maxy: usize = @intFromFloat(@round(@max(0, @min(fcanvas_height, @reduce(.Max, @Vector(4, f32){ p0y, p1y, p2y, p3y })))));
    // raw pointers to avoid slice overhead in the loop
    const dest_ptr: [*]u32 = self.buffer.ptr;
    const src_ptr: [*]const u32 = @ptrCast(@alignCast(image.data.ptr));
    // Go through the bounding box and transform to image space
    for (miny..maxy) |j| {
      for (minx..maxx) |i| {
        const ux, const uy = applyMatrixToPoint(inverse, @floatFromInt(i), @floatFromInt(j));
        // user → image
        const srcX, const srcY = if (comptime no_stretching) lbl: {
          const tx = (ux - fdx);
          const ty = (uy - fdy);
          const srcX = fsx + tx;
          const srcY = fsy + ty;
          break :lbl .{ srcX, srcY };
        } else lbl: {
          const tx = (ux - fdx) / fdWidth;
          const ty = (uy - fdy) / fdHeight;
          const srcX = fsx + tx * fsWidth;
          const srcY = fsy + ty * fsHeight;
          break :lbl .{ srcX, srcY };
        };
        if (srcX < fsx or srcX >= fsx + fsWidth or
            srcY < fsy or srcY >= fsy + fsHeight) {
          continue;
        }
        // Here we choose the interpolation technique
        const fg_packed = lbl: switch (comptime interpolation_type) {
          .nearest => {
            const ix: usize = @intFromFloat(@round(srcX));
            const iy: usize = @intFromFloat(@round(srcY));
            const image_index = iy * image.width + ix;
            break :lbl src_ptr[image_index];
          },
          .linear => {
            const rx = std.math.modf(srcX);
            const ry = std.math.modf(srcY);
            const floor_ix: u16 = @intFromFloat(rx.ipart);
            const floor_iy: u16 = @intFromFloat(ry.ipart);
            const ceil_ix = floor_ix + 1;
            const ceil_iy = floor_iy + 1;
            // Operate on RGBA value in SIMD vectors
            var c00: @Vector(4, f32) = @floatFromInt(@as(@Vector(4, u8), @bitCast(src_ptr[floor_ix + floor_iy * image.width])));
            var c01: @Vector(4, f32) = @floatFromInt(@as(@Vector(4, u8), @bitCast(src_ptr[ceil_ix + floor_iy * image.width])));
            var c10: @Vector(4, f32) = @floatFromInt(@as(@Vector(4, u8), @bitCast(src_ptr[floor_ix + ceil_iy * image.width])));
            var c11: @Vector(4, f32) = @floatFromInt(@as(@Vector(4, u8), @bitCast(src_ptr[ceil_ix + ceil_iy * image.width])));
            if (comptime premultiplyAlpha) {
              // Pre-multiply alpha
              c00 = premulAlpha(c00, c00[3] / 255.0);
              c01 = premulAlpha(c01, c01[3] / 255.0);
              c10 = premulAlpha(c10, c10[3] / 255.0);
              c11 = premulAlpha(c11, c11[3] / 255.0);
            }
            const c0: @Vector(4, f32) = c00 * @as(@Vector(4, f32), @splat(1 - rx.fpart)) + c01 * @as(@Vector(4, f32), @splat(rx.fpart));
            const c1: @Vector(4, f32) = c10 * @as(@Vector(4, f32), @splat(1 - rx.fpart)) + c11 * @as(@Vector(4, f32), @splat(rx.fpart));
            var c: @Vector(4, f32) = c0 * @as(@Vector(4, f32), @splat(1 - ry.fpart)) + c1 * @as(@Vector(4, f32), @splat(ry.fpart));
            if (comptime premultiplyAlpha) {
              if (c[3] > 1.0) {
                // Un-premultiply
                c = premulAlpha(c, 255.0 / c[3]);
              } else {
                c = @splat(0);
              }
            }
            const res: @Vector(4, u8) = @intFromFloat(c);
            break :lbl @as(u32, @bitCast(res));
          },
          // Made with Gemini Pro 3
          .cubic => {
            const rx = std.math.modf(srcX);
            const ry = std.math.modf(srcY);

            // Calculate base integer coordinates
            const base_x: isize = @intFromFloat(rx.ipart);
            const base_y: isize = @intFromFloat(ry.ipart);
            const w: isize = @intCast(image.width);
            const h: isize = @intCast(image.height);

            // Calculate Catmull-Rom weights for X
            // Weights are calculated based on fractional part (fx) for points: -1, 0, 1, 2
            const fx = rx.fpart;
            const fx2 = fx * fx;
            const fx3 = fx2 * fx;
            const wx0: @Vector(4, f32) = @splat(-0.5 * fx3 + fx2 - 0.5 * fx);
            const wx1: @Vector(4, f32) = @splat( 1.5 * fx3 - 2.5 * fx2 + 1.0);
            const wx2: @Vector(4, f32) = @splat(-1.5 * fx3 + 2.0 * fx2 + 0.5 * fx);
            const wx3: @Vector(4, f32) = @splat( 0.5 * fx3 - 0.5 * fx2);

            // Calculate Catmull-Rom weights for Y
            const fy = ry.fpart;
            const fy2 = fy * fy;
            const fy3 = fy2 * fy;
            const wy = [4]@Vector(4, f32){
                @splat(-0.5 * fy3 + fy2 - 0.5 * fy),
                @splat( 1.5 * fy3 - 2.5 * fy2 + 1.0),
                @splat(-1.5 * fy3 + 2.0 * fy2 + 0.5 * fy),
                @splat( 0.5 * fy3 - 0.5 * fy2),
            };

            // Pre-calculate clamped X indices to handle image edges safely
            // (Accessing index -1 or width+1 would crash or read garbage)
            const ix0 = @as(usize, @intCast(std.math.clamp(base_x - 1, 0, w - 1)));
            const ix1 = @as(usize, @intCast(std.math.clamp(base_x,     0, w - 1)));
            const ix2 = @as(usize, @intCast(std.math.clamp(base_x + 1, 0, w - 1)));
            const ix3 = @as(usize, @intCast(std.math.clamp(base_x + 2, 0, w - 1)));

            var color_acc: @Vector(4, f32) = @splat(0);

            // Convolve separate rows (Y axis loop)
            inline for (0..4) |k| {
                // Clamp Y index
                const iy = std.math.clamp(base_y - 1 + @as(isize, k), 0, h - 1);
                const row_offset = @as(usize, @intCast(iy)) * image.width;

                // Fetch row pixels
                var c0: @Vector(4, f32) = @floatFromInt(@as(@Vector(4, u8), @bitCast(src_ptr[ix0 + row_offset])));
                var c1: @Vector(4, f32) = @floatFromInt(@as(@Vector(4, u8), @bitCast(src_ptr[ix1 + row_offset])));
                var c2: @Vector(4, f32) = @floatFromInt(@as(@Vector(4, u8), @bitCast(src_ptr[ix2 + row_offset])));
                var c3: @Vector(4, f32) = @floatFromInt(@as(@Vector(4, u8), @bitCast(src_ptr[ix3 + row_offset])));
                if (comptime premultiplyAlpha) {
                  // Pre-multiply alpha
                  c0 = premulAlpha(c0, c0[3] / 255.0);
                  c1 = premulAlpha(c1, c1[3] / 255.0);
                  c2 = premulAlpha(c2, c2[3] / 255.0);
                  c3 = premulAlpha(c3, c3[3] / 255.0);
                }

                // Interpolate horizontally (X axis)
                const row_val = c0 * wx0 + c1 * wx1 + c2 * wx2 + c3 * wx3;

                // Accumulate vertically (Y axis)
                color_acc += row_val * wy[k];
            }

            if (comptime premultiplyAlpha) {
              if (color_acc[3] > 1.0) {
                // Un-premultiply
                color_acc = premulAlpha(color_acc, 255.0 / color_acc[3]);
              } else {
                color_acc = @splat(0);
              }
            }

            // Clamp results to [0, 255] because cubic interpolation can overshoot
            const clamped = std.math.clamp(color_acc, @as(@Vector(4, f32), @splat(0)), @as(@Vector(4, f32), @splat(255)));

            const c: @Vector(4, u8) = @intFromFloat(clamped);
            break :lbl @as(u32, @bitCast(c));
          },
          // else => @panic(std.fmt.comptimePrint("interpolation type {} not implemented", .{ interpolation_type })),
        };

        const buffer_index = j * self.width + i;
        if (comptime alpha) {
          const bg_packed = dest_ptr[buffer_index];
          const result_8 = blend(bg_packed, fg_packed);
          // Store full pixel (Single instruction store)
          dest_ptr[buffer_index] = result_8;
        } else {
          dest_ptr[buffer_index] = fg_packed;
        }
      }
    }
  }

  fn premulAlpha(v: @Vector(4, f32), a: f32) @Vector(4, f32) {
      return v * @Vector(4, f32){ a, a, a, 1.0 };
  }

  fn blend(bg_packed: u32, fg_packed: u32) u32 {
    // We expand to u16 immediately to prevent overflow during multiply
    const fg: @Vector(4, u16) = @intCast(@as(@Vector(4, u8), @bitCast(fg_packed)));
    const bg: @Vector(4, u16) = @intCast(@as(@Vector(4, u8), @bitCast(bg_packed)));
    const alpha_val = fg[3];
    // Optimization: If alpha is 0, skip. If 255, simple copy.
    if (alpha_val == 0) return bg_packed;
    if (alpha_val == 255) {
        return fg_packed;
    }

    const a: @Vector(4, u16) = @splat(alpha_val);
    const max: @Vector(4, u16) = @splat(255);
    const inv_a = max - a;
    // (bg * inv_a + fg * a) / 255
    const tmp = bg * inv_a + fg * a;
    const result_16 = (tmp + @as(@Vector(4, u16), @splat(1)) + (tmp >> @as(@Vector(4, u8), @splat(8)))) >> @as(@Vector(4, u8), @splat(8));

    const result_8: @Vector(4, u8) = @intCast(result_16);
    return @bitCast(result_8);
  }

  // SIMD version of innerDrawImage3. 25% faster than the regular version
  // Made with Gemini Pro 3
  // pub fn innerDrawImage3(self: *Self, comptime no_stretching: bool, comptime alpha: bool,
  //   comptime interpolation_type: InterpolationType,
  //   image: ImageData, sx: i32, sy: i32, sWidth: usize, sHeight: usize,
  //   dx: i32, dy: i32, dWidth: usize, dHeight: usize) void {
  //   const fsx: f32 = @floatFromInt(sx);
  //   const fsy: f32 = @floatFromInt(sy);
  //   const fsWidth: f32 = @floatFromInt(sWidth);
  //   const fsHeight: f32 = @floatFromInt(sHeight);
  //   const fdx: f32 = @floatFromInt(dx);
  //   const fdy: f32 = @floatFromInt(dy);
  //   const fdWidth: f32 = @floatFromInt(dWidth);
  //   const fdHeight: f32 = @floatFromInt(dHeight);
  //   const fcanvas_width: f32 = @floatFromInt(self.width);
  //   const fcanvas_height: f32 = @floatFromInt(self.height);
  //   _ = interpolation_type; // TODO
  //   // Inverse the matrix transform
  //   var inverse: [6]f32 = .{ 0 } ** 6;
  //   self.getInverse(&inverse) catch @panic("non invertible transform");
  //   // Compute the destination parallelogram
  //   const p0x, const p0y = self.transformPoint(fdx, fdy);
  //   const p1x, const p1y = self.transformPoint(fdx + fdWidth, fdy);
  //   const p2x, const p2y = self.transformPoint(fdx, fdy + fdHeight);
  //   const p3x, const p3y = self.transformPoint(fdx + fdWidth, fdy + fdHeight);
  //   // Compute the bounding box taking into account the canvas dimentions
  //   const minx: usize = @intFromFloat(@round(@max(0, @min(fcanvas_width, @reduce(.Min, @Vector(4, f32){ p0x, p1x, p2x, p3x })))));
  //   const maxx: usize = @intFromFloat(@round(@max(0, @min(fcanvas_width, @reduce(.Max, @Vector(4, f32){ p0x, p1x, p2x, p3x })))));
  //   const miny: usize = @intFromFloat(@round(@max(0, @min(fcanvas_height, @reduce(.Min, @Vector(4, f32){ p0y, p1y, p2y, p3y })))));
  //   const maxy: usize = @intFromFloat(@round(@max(0, @min(fcanvas_height, @reduce(.Max, @Vector(4, f32){ p0y, p1y, p2y, p3y })))));
  //   // raw pointers to avoid slice overhead in the loop
  //   const dest_ptr: [*]u32 = self.buffer.ptr;
  //   const src_ptr: [*]const u32 = @ptrCast(@alignCast(image.data.ptr));

  //   // --- SIMD PRE-CALCULATION START ---
  //   // Prepare constants to replace 'applyMatrixToPoint' and 'no_stretching' logic inside the loop
  //   const scale_x = if (comptime no_stretching) 1.0 else fsWidth / fdWidth;
  //   const scale_y = if (comptime no_stretching) 1.0 else fsHeight / fdHeight;

  //   // Combined coefficients: src = (Inverse * Screen - dest_offset) * scale + src_offset
  //   // Assuming standard 2D matrix layout: x' = inv[0]x + inv[2]y + inv[4]
  //   const A = inverse[0] * scale_x;
  //   const B = inverse[2] * scale_x;
  //   const C = (inverse[4] - fdx) * scale_x + fsx;

  //   const D = inverse[1] * scale_y;
  //   const E = inverse[3] * scale_y;
  //   const F = (inverse[5] - fdy) * scale_y + fsy;

  //   const Vec4f = @Vector(4, f32);
  //   const Vec4i = @Vector(4, i32);
  //   const vec_offset = Vec4f{ 0, 1, 2, 3 };
  //   const vec_A = @as(Vec4f, @splat(A));
  //   const vec_D = @as(Vec4f, @splat(D));
  //   // --- SIMD PRE-CALCULATION END ---

  //   // Go through the bounding box and transform to image space
  //   for (miny..maxy) |j| {
  //     // Pre-calculate row base components
  //     const fj = @as(f32, @floatFromInt(j));
  //     const row_base_x = B * fj + C;
  //     const row_base_y = E * fj + F;

  //     var i = minx;
  //     // Vectorized Loop (4 pixels at a time)
  //     while (i + 4 <= maxx) : (i += 4) {
  //       const fi = @as(f32, @floatFromInt(i));

  //       // Calculate Screen Coordinates: {i, i+1, i+2, i+3}
  //       const vec_scr_x = @as(Vec4f, @splat(fi)) + vec_offset;

  //       // Apply Matrix and Scaling (SIMD)
  //       const vec_src_x = vec_A * vec_scr_x + @as(Vec4f, @splat(row_base_x));
  //       const vec_src_y = vec_D * vec_scr_x + @as(Vec4f, @splat(row_base_y));

  //       // Convert to Integer (SIMD) - Rounding as per original code
  //       const vec_ix = @as(Vec4i, @intFromFloat(@round(vec_src_x)));
  //       const vec_iy = @as(Vec4i, @intFromFloat(@round(vec_src_y)));

  //       // Extract to array to handle memory access
  //       const arr_ix: [4]i32 = vec_ix;
  //       const arr_iy: [4]i32 = vec_iy;
  //       const arr_src_x: [4]f32 = vec_src_x;
  //       const arr_src_y: [4]f32 = vec_src_y;

  //       inline for (0..4) |k| {
  //         const srcX = arr_src_x[k];
  //         const srcY = arr_src_y[k];

  //         if (srcX >= fsx and srcX < fsx + fsWidth and
  //           srcY >= fsy and srcY < fsy + fsHeight) {

  //           const ix: usize = @intCast(arr_ix[k]);
  //           const iy: usize = @intCast(arr_iy[k]);

  //           const buffer_index = j * self.width + (i + k);
  //           const image_index = iy * image.width + ix;
  //           const fg_packed = src_ptr[image_index];

  //           if (comptime alpha) {
  //             const bg_packed = dest_ptr[buffer_index];
  //             const result_8 = blend(bg_packed, fg_packed);
  //             dest_ptr[buffer_index] = result_8;
  //           } else {
  //             dest_ptr[buffer_index] = fg_packed;
  //           }
  //         }
  //       }
  //     }

  //     // Handle remaining pixels (Scalar fallback)
  //     while (i < maxx) : (i += 1) {
  //       const fi = @as(f32, @floatFromInt(i));
  //       const srcX = A * fi + row_base_x;
  //       const srcY = D * fi + row_base_y;

  //       if (srcX >= fsx and srcX < fsx + fsWidth and
  //         srcY >= fsy and srcY < fsy + fsHeight) {

  //         const ix: usize = @intFromFloat(@round(srcX));
  //         const iy: usize = @intFromFloat(@round(srcY));
  //         const buffer_index = j * self.width + i;
  //         const image_index = iy * image.width + ix;
  //         const fg_packed = src_ptr[image_index];

  //         if (comptime alpha) {
  //           const bg_packed = dest_ptr[buffer_index];
  //           const result_8 = blend(bg_packed, fg_packed);
  //           dest_ptr[buffer_index] = result_8;
  //         } else {
  //           dest_ptr[buffer_index] = fg_packed;
  //         }
  //       }
  //     }
  //   }
  // }

  //
  // Private functions
  //

  fn getInterpolationType(self: Self) InterpolationType {
    if (self.imageSmoothingEnabled) {
      if (std.mem.eql(u8, self.imageSmoothingQuality, "low")) {
        return InterpolationType.nearest;
      } else if (std.mem.eql(u8, self.imageSmoothingQuality, "high")) {
        return InterpolationType.cubic;
      }
      return InterpolationType.linear;
    }
    return InterpolationType.nearest;
  }

  /// Draws a point without transformation. To be used by private functions.
  inline fn internalPlot(self: Self, x: i16, y: i16, acolor: u32) void {
    if (x >= 0 and x < self.width and y >= 0 and y < self.height) {
      self.buffer[@as(u32, @intCast(y)) * self.width + @as(u32, @intCast(x))] = acolor;
    }
  }

  /// Fills a polygon defined by `vertices` using the scan-line algorithm.
  fn fillPolygon(self: Self, allocator: std.mem.Allocator, vertices: []const [2]f32, afillStyle: u32) void {
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
      fn initializeEdges(vrtces: []const [2]f32, edges: *std.ArrayList(Edge)) void {
        for (vrtces, 0..) |p1, i| {
          const p2 = vrtces[(i + 1) % vrtces.len];

          // Round y-coordinates to determine integer scan-lines.
          const y0_rounded = @round(p1[1]);
          const y1_rounded = @round(p2[1]);

          const min_y: i16 = @intFromFloat(@min(y0_rounded, y1_rounded));
          const max_y: i16 = @intFromFloat(@max(y0_rounded, y1_rounded));
          const x_at_min_y = if (p1[1] <= p2[1]) p1[0] else p2[0];
          const invslope = calculateInvSlope(p1[0], p1[1], p2[0], p2[1]);

          edges.append(Edge{
            .min_y = min_y,
            .max_y = max_y,
            .x_at_min_y = x_at_min_y,
            .invslope = invslope,
          }) catch @panic("OOM");
        }
      }

      /// Define a function to initialize the global edge table.
      fn initializeGlobalEdgeTable(edges: []const Edge, global_edge_table: *std.ArrayList(Edge)) void {
        global_edge_table.appendSlice(edges) catch @panic("OOM");

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
        active_edge_table: *std.ArrayList(Edge),
      ) void {
        active_edge_table.clearRetainingCapacity();
        for (source_table) |edge| {
          if (edge.min_y <= scan_line and edge.max_y > scan_line) {
            active_edge_table.append(edge) catch @panic("OOM");
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

    var edge_list = std.ArrayList(Fns.Edge).init(allocator);
    defer edge_list.deinit();
    Fns.initializeEdges(vertices, &edge_list);

    var global_edge_table = std.ArrayList(Fns.Edge).init(allocator);
    defer global_edge_table.deinit();
    Fns.initializeGlobalEdgeTable(edge_list.items, &global_edge_table);

    if (global_edge_table.items.len == 0) return;

    // Initialize scan-line and active edge table.
    var scan_line: i16 = global_edge_table.items[0].min_y;
    var active_edge_table = std.ArrayList(Fns.Edge).init(allocator);
    defer active_edge_table.deinit();
    Fns.getActiveEdgeTable(global_edge_table.items, scan_line, &active_edge_table);

    // Iterate over each scan-line until active edge table is empty.
    while (active_edge_table.items.len > 0) {
      std.log.debug("scanline {}", .{ active_edge_table.items.len });
      // Draw pixels between x-values of odd and even parity edge pairs.
      var i: usize = 0;
      while (i + 1 < active_edge_table.items.len) : (i += 2) {
        const even_edge = active_edge_table.items[i];
        const odd_edge = active_edge_table.items[i + 1];
        std.log.debug("even_edge {}, odd_edge {}", .{ even_edge, odd_edge });

        const x_start_f = asf32(scan_line - even_edge.min_y) * even_edge.invslope + even_edge.x_at_min_y;
        const x_end_f = asf32(scan_line - odd_edge.min_y) * odd_edge.invslope + odd_edge.x_at_min_y;

        var x_start = @round(x_start_f);
        const x_end = @round(x_end_f);

        std.log.debug("x_start {} x_end {}", .{ x_start, x_end });
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
    return applyMatrixToPoint(self.transformMatrix, x, y);
  }

  fn applyMatrixToPoint(m: [6]f32, x: f32, y: f32) struct { f32, f32 } {
    return .{
      m[_a] * x + m[_c] * y + m[_e],
      m[_b] * x + m[_d] * y + m[_f],
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
    self.internalPlot(startx, starty, self.strokeStyle);
    if (tDeltaX > tDeltaY) {
      // start value represents a half step in Y direction
      tError = tDeltaYTimes2 - tDeltaX;
      while (startx != endx) {
        // step in main direction
        startx += tStepX;
        if (tError >= 0) {
          if (aOverlap & @intFromEnum(Overlap.LINE_OVERLAP_MAJOR) != 0) {
            // draw pixel in main direction before changing
            self.internalPlot(startx, starty, self.strokeStyle);
          }
          // change Y
          starty += tStepY;
          if (aOverlap & @intFromEnum(Overlap.LINE_OVERLAP_MINOR) != 0) {
            // draw pixel in minor direction before changing
            self.internalPlot(startx - tStepX, starty, self.strokeStyle);
          }
          tError -= tDeltaXTimes2;
        }
        tError += tDeltaYTimes2;
        self.internalPlot(startx, starty, self.strokeStyle);
      }
    } else {
      tError = tDeltaXTimes2 - tDeltaY;
      while (starty != endy) {
        starty += tStepY;
        if (tError >= 0) {
          if (aOverlap & @intFromEnum(Overlap.LINE_OVERLAP_MAJOR) != 0) {
            // draw pixel in main direction before changing
            self.internalPlot(startx, starty, self.strokeStyle);
          }
          startx += tStepX;
          if (aOverlap & @intFromEnum(Overlap.LINE_OVERLAP_MINOR) != 0) {
            // draw pixel in minor direction before changing
            self.internalPlot(startx, starty - tStepY, self.strokeStyle);
          }
          tError -= tDeltaYTimes2;
        }
        tError += tDeltaXTimes2;
        self.internalPlot(startx, starty, self.strokeStyle);
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

  pub fn lerp(start: @Vector(2, f32), end: @Vector(2, f32), t: f32) @Vector(2, f32) {
    return start + (end - start) * @as(@Vector(2, f32), @splat(t));
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
  // {
  //   var vertices = [_][2]f32{ [2]f32{ 1.0, 1.0 }, [2]f32{ 3.0, 1.0 }, [2]f32{ 3.0, 3.0 }, [2]f32{ 1.0, 3.0 } };
  //   var ctx = try DrawContext.init(allocator, 5, 5);
  //   defer ctx.deinit(allocator);
  //   ctx.fillPolygon(allocator, &vertices, 1);
  //   // try printBuffer(ctx);
  //   // FIXME: We shouldn't paint the edge of the polygon
  //   try std.testing.expectEqualSlices(u32, &.{
  //     // 0, 0, 0, 0, 0,
  //     // 0, 0, 0, 0, 0,
  //     // 0, 0, 1, 0, 0,
  //     // 0, 0, 0, 0, 0,
  //     // 0, 0, 0, 0, 0,
  //     0, 0, 0, 0, 0,
  //     0, 1, 1, 0, 0,
  //     0, 1, 1, 0, 0,
  //     0, 0, 0, 0, 0,
  //     0, 0, 0, 0, 0,
  //   }, ctx.buffer);
  // }
  {
    var vertices = [_][2]f32{ .{ 4, 0 }, .{ 1, 3 }, .{ 3, 3 }, .{ 3, 2 } };
    var ctx = try DrawContext.init(allocator, 5, 5);
    defer ctx.deinit(allocator);
    ctx.fillPolygon(allocator, &vertices, 1);
    try printBuffer(ctx);
    try std.testing.expectEqualSlices(u32, &.{
      // 0, 0, 0, 0, 0,
      // 0, 0, 0, 0, 0,
      // 0, 0, 1, 0, 0,
      // 0, 0, 0, 0, 0,
      // 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0,
      0, 0, 0, 1, 0,
      0, 0, 1, 1, 0,
      0, 1, 1, 1, 0,
      0, 0, 0, 0, 0,
    }, ctx.buffer);
  }
}

test "fillCircle" {
  std.testing.log_level = .debug;
  const allocator = std.testing.allocator;

  var ctx = try DrawContext.init(allocator, 5, 5);
  defer ctx.deinit(allocator);
  ctx.fillStyle = 1;
  ctx.fillCircle(2, 2, 2);
  // try printBuffer(ctx);
  try std.testing.expectEqualSlices(u32, &.{
    0, 0, 1, 0, 0,
    0, 1, 1, 1, 0,
    1, 1, 1, 1, 1,
    0, 1, 1, 1, 0,
    0, 0, 1, 0, 0,
  }, ctx.buffer);
}

test "getInverse" {
  std.testing.log_level = .debug;
  const allocator = std.testing.allocator;

  var ctx = try DrawContext.init(allocator, 5, 5);
  defer ctx.deinit(allocator);

  var inverse: [6]f32 = .{ 0 } ** 6;
  try ctx.getInverse(&inverse);
  try std.testing.expectEqual(.{ 1, 0, 0, 1, 0, 0 }, inverse);

  ctx.translate(2, 0);
  try ctx.getInverse(&inverse);
  try std.testing.expectEqual(.{ 1, 0, 0, 1, -2, 0 }, inverse);

  ctx.rotate(std.math.pi);
  try ctx.getInverse(&inverse);
  try std.testing.expect(-1 - inverse[0] < 0.000001);
  try std.testing.expect( 0 - inverse[1] < 0.000001);
  try std.testing.expect( 0 - inverse[2] < 0.000001);
  try std.testing.expect(-1 - inverse[3] < 0.000001);
  try std.testing.expect( 2 - inverse[4] < 0.000001);
  try std.testing.expect( 0 - inverse[5] < 0.000001);
}

test "zoom" {
  std.testing.log_level = .debug;
  const allocator = std.testing.allocator;

  var ctx = try DrawContext.init(allocator, 5, 5);
  defer ctx.deinit(allocator);
  ctx.translate(2, 2);
  ctx.scale(2, 2);
  ctx.translate(-2, -2);
  try std.testing.expectEqual(.{ 2e0, 0e0, 0e0, 2e0, -2e0, -2e0 }, ctx.getTransform());
}

test "line" {
  std.testing.log_level = .debug;
  const allocator = std.testing.allocator;

  var ctx = try DrawContext.init(allocator, 5, 5);
  defer ctx.deinit(allocator);
  ctx.scale(2, 2);
  ctx.strokeStyle = 1;
  ctx.line(0, 0, 5, 5);
  try printBuffer(ctx);
  try std.testing.expectEqualSlices(u32, &.{
    1, 0, 0, 0, 0,
    0, 1, 0, 0, 0,
    0, 0, 1, 0, 0,
    0, 0, 0, 1, 0,
    0, 0, 0, 0, 1,
  }, ctx.buffer);
  // try std.testing.expectEqual(.{ 2e0, 0e0, 0e0, 2e0, -2e0, -2e0 }, ctx.getTransform());
}
