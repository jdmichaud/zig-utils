const std = @import("std");

const fontfile = @embedFile("dos_8x8_font_white.pbm");

pub fn DrawContext(comptime pwidth: u32, comptime pheight: u32) type {
    return struct {
        const Self = @This();

        pub const contextWidth = pwidth;
        pub const contextHeight = pheight;

        pub var buffer = [_]u32{0} ** (pwidth * pheight);

        pub var fillStyle: u32 = 0;
        pub var strokeStyle: u32 = 0;
        pub var lineWidth: u32 = 0;

        var _transform: [6]f32 = .{ 1, 0, 0, 1, 0, 0 };
        const _a = 0;
        const _b = 1;
        const _c = 2;
        const _d = 3;
        const _e = 4;
        const _f = 5;

        var stack: [6]f32 = .{ 1, 0, 0, 1, 0, 0 };

        // pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
        //   allocator.free(self.font);
        // }

        // resets (overrides) the current transformation to the identity matrix, and
        // then invokes a transformation described by the arguments of this method.
        // This lets you scale, rotate, translate (move), and skew the context.
        //                                              a c e
        // The transformation matrix is described by: [ b d f ]
        //                                              0 0 1
        pub fn setTransform(a: f32, b: f32, c: f32, d: f32, e: f32, f: f32) void {
            _transform = .{ a, b, c, d, e, f };
        }
        // multiplies the current transformation with the matrix described by the
        // arguments of this method. This lets you scale, rotate, translate (move),
        // and skew the context.
        pub fn transform(a: f32, b: f32, c: f32, d: f32, e: f32, f: f32) void {
            _transform[_a] = _transform[_a] * a + _transform[_c] * b;
            _transform[_b] = _transform[_b] * a + _transform[_d] * b;
            _transform[_c] = _transform[_a] * c + _transform[_c] * d;
            _transform[_d] = _transform[_b] * c + _transform[_d] * d;
            _transform[_e] = _transform[_a] * e + _transform[_c] * f + _transform[_e] * 1;
            _transform[_f] = _transform[_b] * e + _transform[_d] * f + _transform[_f] * 1;
        }
        // retrieves the current transformation matrix being applied to the context.
        pub fn getTransform() [6]f32 {
            return _transform;
        }
        // Saves the entire state of the canvas by pushing the current state onto a
        // stack.
        // ⚠️ Only one level of stack for now.
        pub fn save() !void {
            stack = _transform;
        }
        // Restores the most recently saved canvas state by popping the top entry
        // in the drawing state stack. If there is no saved state, this method does
        // nothing.
        pub fn restore() void {
            _transform = stack;
        }
        // resets the rendering context to its default state, allowing it to be
        // reused for drawing something else without having to explicitly reset all
        // the properties.
        pub fn reset() void {
            _transform = .{ 1, 0, 0, 1, 0, 0 };
        }
        // adds a translation transformation to the current matrix.
        pub fn translate(x: f32, y: f32) void {
            _transform[_e] += x;
            _transform[_f] += y;
        }
        // fill a rectangle with the fillStyle
        pub fn fillRect(x: i16, y: i16, width: i16, height: i16) void {
            const xc, const yc = transformPoint(x, y);
            const widthc, const heightc = transformVector(width, height);
            fillPolygon(&.{
                &.{ xc, yc },
                &.{ xc + widthc, yc },
                &.{ xc + widthc, yc + heightc },
                &.{ xc, yc + heightc },
            }, fillStyle);
        }
        // erases the pixels in a rectangular area by setting them to transparent
        // black.
        pub fn clearRect(x: i16, y: i16, width: i16, height: i16) void {
            const previousFillStyle = fillStyle;
            fillStyle = 0xFF000000;
            fillRect(x, y, width, height);
            fillStyle = previousFillStyle;
        }
        // Renders a rectangle with a starting point is at (x, y) and whose size is
        // specified by width and height.
        pub fn strokeRect(x: i16, y: i16, width: i16, height: i16) void {
            line(x, y, x + width, y);
            line(x + width, y, x + width, y + height);
            line(x + width, y + height, x, y + height);
            line(x, y + height, x, y);
        }
        // Draws a line.
        pub fn line(startx: i16, starty: i16, endx: i16, endy: i16) void {
            // std.debug.assert(startx >= 0 and starty >= 0 and endx >= 0 and endy >= 0);

            // This is commented because, although debug performances are much better,
            // release performances are worst!
            const builtin = @import("builtin");
            if (builtin.mode == .Debug) {
                // This whole block is an optimization for vertical and horizontal line
                const ux =
                    _transform[_a] * @as(f32, @floatFromInt(startx)) +
                    _transform[_c] * @as(f32, @floatFromInt(starty)) + _transform[_e];
                var uy =
                    _transform[_b] * @as(f32, @floatFromInt(startx)) +
                    _transform[_d] * @as(f32, @floatFromInt(starty)) + _transform[_f];
                const vx =
                    _transform[_a] * @as(f32, @floatFromInt(endx)) +
                    _transform[_c] * @as(f32, @floatFromInt(endy)) + _transform[_e];
                var vy =
                    _transform[_b] * @as(f32, @floatFromInt(endx)) +
                    _transform[_d] * @as(f32, @floatFromInt(endy)) + _transform[_f];

                if (ux > 0 and uy > 0 and vx > 0 and vy > 0 and
                    ux < contextWidth and uy < contextHeight and vx < contextWidth and vy < contextHeight)
                {
                    // If a line is entirely in the canvas
                    if (ux == vx) {
                        // vertical line
                        const x: u16 = @intFromFloat(ux);
                        if (uy > vy) std.mem.swap(@TypeOf(uy), &uy, &vy);
                        var y: u16 = @intFromFloat(uy);
                        while (y < @as(u16, @intFromFloat(vy))) : (y += 1) {
                            buffer[y * contextWidth + x] = strokeStyle;
                        }
                        return;
                    } else if (uy == vy) {
                        // horizontal line
                        var startBuffer = @as(u16, @intFromFloat(uy)) * contextWidth + @as(u16, @intFromFloat(ux));
                        var endBuffer = @as(u16, @intFromFloat(vy)) * contextWidth + @as(u16, @intFromFloat(vx));
                        if (startBuffer > endBuffer) std.mem.swap(@TypeOf(startBuffer), &startBuffer, &endBuffer);
                        @memset(buffer[startBuffer..endBuffer], strokeStyle);
                        return;
                    }
                }
            }
            // Otherwise, we use a general but slow algorithm
            drawThickLine(startx, starty, endx, endy);
            // drawLineOverlap(startx, starty, endx, endy, 0);
            // drawLineWu(startx, starty, endx, endy, 0);
        }
        // Draws a point.
        pub inline fn plot(x: i16, y: i16, acolor: u32) void {
            // std.log.debug("plot x {} y {} width {} height {} index {} buffer.len {}", .{
            //   x, y, contextWidth, contextHeight,
            //   @as(u16, @bitCast(y)) * width + @as(u16, @bitCast(x)),
            //   buffer.len,
            // });
            const vx = _transform[_a] * @as(f32, @floatFromInt(x)) + _transform[_c] * @as(f32, @floatFromInt(y)) + _transform[_e];
            const vy = _transform[_b] * @as(f32, @floatFromInt(x)) + _transform[_d] * @as(f32, @floatFromInt(y)) + _transform[_f];
            if (vx >= 0 and vx < contextWidth and vy >= 0 and vy < contextHeight) {
                buffer[@as(u16, @intFromFloat(vy)) * contextWidth + @as(u16, @intFromFloat(vx))] = acolor;
            }
        }
        inline fn plot_internal(x: i16, y: i16, acolor: u32) void {
            buffer[y * contextWidth + x] = acolor;
        }
        // Writes text at the specified position. x and y specifies the top left
        // corner of the text box to be printed.
        pub fn fillText(text: []const u8, x: i16, y: i16) void {
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
                            plot(x + @as(i16, @intCast(i)) + @as(i16, @intCast(cindex * 8)), y + @as(i16, @intCast(j)), fillStyle);
                        }
                    }
                }
            }
        }

        // Private functions

        /// Transform a point according to the currently set transform matrix in the context.
        fn transformPoint(x: f32, y: f32) struct { f32, f32 } {
            return .{
                _transform[_a] * x + _transform[_c] * y + _transform[_e],
                _transform[_b] * x + _transform[_d] * y + _transform[_f],
            };
        }

        /// Transform a vector according to the currently set transform matrix in the context.
        fn transformVector(x: f32, y: f32) struct { f32, f32 } {
            return .{
                _transform[_a] * x + _transform[_c] * y,
                _transform[_b] * x + _transform[_d] * y,
            };
        }

        // Modified Bresenham draw(line) with optional overlap. Required for drawThickLine().
        // Overlap draws additional pixel when changing minor direction. For standard bresenham overlap, choose LINE_OVERLAP_NONE (0).
        //
        //  Sample line:
        //
        //    00+
        //     -0000+
        //         -0000+
        //             -00
        //
        //  0 pixels are drawn for normal line without any overlap LINE_OVERLAP_NONE
        //  + pixels are drawn if LINE_OVERLAP_MAJOR
        //  - pixels are drawn if LINE_OVERLAP_MINOR

        const Overlap = enum(u8) {
            LINE_OVERLAP_NONE = 0,
            LINE_OVERLAP_MINOR = 1,
            LINE_OVERLAP_MAJOR = 2,
            LINE_OVERLAP_BOTH = 3,
        };

        fn drawLineOverlap(pstartx: i16, pstarty: i16, endx: i16, endy: i16, aOverlap: u8) void {
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
            plot(startx, starty, strokeStyle);
            if (tDeltaX > tDeltaY) {
                // start value represents a half step in Y direction
                tError = tDeltaYTimes2 - tDeltaX;
                while (startx != endx) {
                    // step in main direction
                    startx += tStepX;
                    if (tError >= 0) {
                        if (aOverlap & @intFromEnum(Overlap.LINE_OVERLAP_MAJOR) != 0) {
                            // draw pixel in main direction before changing
                            plot(startx, starty, strokeStyle);
                        }
                        // change Y
                        starty += tStepY;
                        if (aOverlap & @intFromEnum(Overlap.LINE_OVERLAP_MINOR) != 0) {
                            // draw pixel in minor direction before changing
                            plot(startx - tStepX, starty, strokeStyle);
                        }
                        tError -= tDeltaXTimes2;
                    }
                    tError += tDeltaYTimes2;
                    plot(startx, starty, strokeStyle);
                }
            } else {
                tError = tDeltaXTimes2 - tDeltaY;
                while (starty != endy) {
                    starty += tStepY;
                    if (tError >= 0) {
                        if (aOverlap & @intFromEnum(Overlap.LINE_OVERLAP_MAJOR) != 0) {
                            // draw pixel in main direction before changing
                            plot(startx, starty, strokeStyle);
                        }
                        startx += tStepX;
                        if (aOverlap & @intFromEnum(Overlap.LINE_OVERLAP_MINOR) != 0) {
                            // draw pixel in minor direction before changing
                            plot(startx, starty - tStepY, strokeStyle);
                        }
                        tError -= tDeltaYTimes2;
                    }
                    tError += tDeltaXTimes2;
                    plot(startx, starty, strokeStyle);
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

        fn shadeColor(R: f32, G: f32, B: f32, ratio: f32) u32 {
            return @as(u32, @intFromFloat(R * ratio)) |
                @as(u32, @intFromFloat(G * ratio)) << 8 |
                @as(u32, @intFromFloat(B * ratio)) << 16 |
                0xFF000000; // Always full alpha
        }

        fn drawLineWu(px0: i16, py0: i16, px1: i16, py1: i16, unused: u8) void {
            std.log.debug("drawLineWu {} {} {} {}", .{ px0, py0, px1, py1 });
            _ = unused;
            var x0: f32 = @floatFromInt(px0);
            var y0: f32 = @floatFromInt(py0);
            var x1: f32 = @floatFromInt(px1);
            var y1: f32 = @floatFromInt(py1);
            const R: f32 = @floatFromInt(strokeStyle & 0x000000FF);
            const G: f32 = @floatFromInt((strokeStyle & 0x0000FF00) >> 8);
            const B: f32 = @floatFromInt((strokeStyle & 0x00FF0000) >> 16);

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
                plot(ypxl1, xpxl1, shadeColor(R, G, B, rfpart(yend) * xgap));
                plot(ypxl1 + 1, xpxl1, shadeColor(R, G, B, fpart(yend) * xgap));
            } else {
                plot(xpxl1, ypxl1, shadeColor(R, G, B, rfpart(yend) * xgap));
                plot(xpxl1, ypxl1 + 1, shadeColor(R, G, B, fpart(yend) * xgap));
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
                plot(ypxl2, xpxl2, shadeColor(R, G, B, rfpart(yend) * xgap));
                plot(ypxl2 + 1, xpxl2, shadeColor(R, G, B, fpart(yend) * xgap));
            } else {
                plot(xpxl2, ypxl2, shadeColor(R, G, B, rfpart(yend) * xgap));
                plot(xpxl2, ypxl2 + 1, shadeColor(R, G, B, fpart(yend) * xgap));
            }

            // main loop
            if (steep) {
                var i = xpxl1 + 1;
                while (i < xpxl2 - 1) {
                    plot(@as(i16, @intFromFloat(@floor(intery))), i, shadeColor(R, G, B, rfpart(intery)));
                    plot(@as(i16, @intFromFloat(@floor(intery))) + 1, i, shadeColor(R, G, B, fpart(intery)));
                    intery = intery + gradient;
                    i += 1;
                }
            } else {
                var i = xpxl1 + 1;
                while (i < xpxl2 - 1) {
                    plot(i, @as(i16, @intFromFloat(@floor(intery))), shadeColor(R, G, B, rfpart(intery)));
                    plot(i, @as(i16, @intFromFloat(@floor(intery))) + 1, shadeColor(R, G, B, fpart(intery)));
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
        fn drawThickLine(pstartx: i16, pstarty: i16, pendx: i16, pendy: i16) void {
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
                    var i = lineWidth / 2;
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
                drawLineOverlap(startx, starty, endx, endy, @intFromEnum(tOverlap));
                // drawLineWu(startx, starty, endx, endy, @intFromEnum(tOverlap));
                // draw lineWidth lines
                tError = tDeltaYTimes2 - tDeltaX;
                var i = lineWidth;
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
                    drawLineOverlap(startx, starty, endx, endy, @intFromEnum(tOverlap));
                    // drawLineWu(startx, starty, endx, endy, @intFromEnum(tOverlap));
                    i -= 1;
                }
            } else {
                // if (we want to draw the original coordinate in the middle of the thick line)
                {
                    tError = tDeltaXTimes2 - tDeltaY;
                    var i = lineWidth / 2;
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
                drawLineOverlap(startx, starty, endx, endy, @intFromEnum(tOverlap));
                // drawLineWu(startx, starty, endx, endy, @intFromEnum(tOverlap));
                tError = tDeltaXTimes2 - tDeltaY;
                var i = lineWidth;
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
                    drawLineOverlap(startx, starty, endx, endy, @intFromEnum(tOverlap));
                    // drawLineWu(startx, starty, endx, endy, @intFromEnum(tOverlap));
                    i -= 1;
                }
            }
        }

        /// Fills a polygon defined by `vertices` using the scan-line algorithm.
        fn fillPolygon(vertices: [2]f32, afillStyle: u32) !void {
            const Fns = struct {
                /// Represents a single edge of the polygon, prepared for the scan-line algorithm.
                /// This version stores 'slope' (dy/dx) to match the original Python code.
                const Edge = struct {
                    min_y: f32,
                    max_y: f32,
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
                fn initializeEdges(vrtces: [2]f32, edges: *std.ArrayListAligned(Edge)) !void {
                    for (vrtces, 0..) |p1, i| {
                        const p2 = vrtces[(i + 1) % vrtces.len];

                        // Round y-coordinates to determine integer scan-lines.
                        const y0_rounded = @round(p1[1]);
                        const y1_rounded = @round(p2[1]);

                        const min_y: i32 = @intFromFloat(@min(y0_rounded, y1_rounded));
                        const max_y: i32 = @intFromFloat(@max(y0_rounded, y1_rounded));
                        const x_at_min_y = if (p1[1] <= p2[1]) p1[0] else p2[0];
                        const invslope = calculateInvSlope(p1[0], p1[1], p2[0], p2[1]);

                        try edges.append(Edge{
                            .min_y = min_y,
                            .max_y = max_y,
                            .x_at_min_y = x_at_min_y,
                            .invslope = invslope,
                        });
                    }
                }

                /// Define a function to initialize the global edge table.
                fn initializeGlobalEdgeTable(edges: []const Edge, global_edge_table: *std.ArrayListAligned(Edge)) !void {
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
                    std.sort.sort(Edge, global_edge_table, {}, sort_context.lessThan);

                    // We eliminate flat edges (slope == 0)
                    var i = 0;
                    var len = global_edge_table.item.len;
                    while (i < len) {
                        if (global_edge_table.item[i].invslope == std.math.inf(f32)) {
                            global_edge_table.orderedRemove(i);
                            len = global_edge_table.item.len;
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
                    active_edge_table: *std.ArrayListAligned(Edge),
                ) !void {
                    for (source_table) |edge| {
                        if (edge.min_y <= scan_line and edge.max_y > scan_line) {
                            try active_edge_table.append(edge);
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
                    std.sort.sort(Edge, active_edge_table.items, {}, sort_context.lessThan);

                    return active_edge_table;
                }
            };

            if (vertices.len < 3) return;

            const edge_buffer: [64]Fns.Edge = undefined; // Preallocate a buffer for edges.
            var edge_list = std.ArrayListAligned.initBuffer(edge_buffer);
            try Fns.initializeEdges(vertices, &edge_list);

            const global_edge_buffer: [64]Fns.Edge = undefined; // Preallocate a buffer for edges.
            var global_edge_table = std.ArrayListAligned.initBuffer(global_edge_buffer);
            try Fns.initializeGlobalEdgeTable(edge_list.items, &global_edge_table);

            if (global_edge_table.len == 0) return;

            // Initialize scan-line and active edge table.
            var scan_line: i32 = global_edge_table[0].min_y;
            var active_edge_table = std.ArrayListAligned.initBuffer(edge_buffer);
            try Fns.getActiveEdgeTable(global_edge_table, scan_line, &active_edge_table);

            // Iterate over each scan-line until active edge table is empty.
            while (active_edge_table.items.len > 0) {
                // Draw pixels between x-values of odd and even parity edge pairs.
                var i: usize = 0;
                while (i + 1 < active_edge_table.items.len) : (i += 2) {
                    const even_edge = active_edge_table.items[i];
                    const odd_edge = active_edge_table.items[i + 1];

                    const scan_line_f: f32 = @floatFromInt(scan_line);

                    const x_start_f = (scan_line_f - even_edge.min_y) * even_edge.invslope + even_edge.x_at_min_y;
                    const x_end_f = (scan_line_f - odd_edge.min_y) * odd_edge.invslope + odd_edge.x_at_min_y;

                    var x_start = @round(x_start_f);
                    const x_end = @round(x_end_f);

                    while (x_start < x_end) {
                        plot(x_start, scan_line, afillStyle);
                        x_start += 1;
                    }
                }

                // Update scan-line and rebuild active edge table for the next line.
                scan_line += 1;
                try Fns.getActiveEdgeTable(global_edge_table, scan_line, &active_edge_table);
            }
        }
    };
}

fn make_polygon(vertices: [][2]f32) void {
    const Fns = struct {
        // compute in which quadrant the vertex belongs to w.r.t. the centroid
        fn quadrant(vertex: [2]f32, centroid: [2]f32) u8 {
            const dx = vertex[0] - centroid[0];
            const dy = vertex[1] - centroid[1];
            if (dx >= 0 and dy >= 0) {
                return 0;
            } else if (dx < 0 and dy >= 0) {
                return 1;
            } else if (dx < 0 and dy < 0) {
                return 2;
            }
            return 3;
        }

        // cross product between two vertices assuming centroid as the origin
        fn cross(v1: [2]f32, v2: [2]f32, centroid: [2]f32) f32 {
            const d1x = v1[0] - centroid[0];
            const d1y = v1[1] - centroid[1];
            const d2x = v2[0] - centroid[0];
            const d2y = v2[1] - centroid[1];
            return d1x * d2y - d1y * d2x;
        }

        // distance from the centroid
        fn dist_sq(vertex: [2]f32, centroid: [2]f32) f32 {
            return std.math.pow(vertex[0] - centroid[0], 2) + std.math.pow(vertex[1] - centroid[1], 2);
        }

        // order by quadrant, then by cross product and then by distance from centroid
        fn cmpfn(context: struct { centroid: [2]f32 }, a: [2]f32, b: [2]f32) f32 {
            const qa = quadrant(a, context.centroid);
            const qb = quadrant(b, context.centroid);
            if (qa != qb) {
                return qa - qb;
            }
            // in same quadrant: sort by cross product
            const c = cross(a, b, context.centroid);
            if (c != 0) {
                return if (c > 0) -1 else 1;
            }
            // Collinear points: sort by distance from centroid
            const d1 = dist_sq(a, context.centroid);
            const d2 = dist_sq(b, context.centroid);
            return d1 - d2;
        }
    };

    var average_x: f32 = 0;
    var average_y: f32 = 0;
    for (vertices) |vertex| {
        average_x += vertex[0];
        average_y += vertex[1];
    }
    average_x /= vertices.len;
    average_y /= vertices.len;
    const centroid: [2]f32 = .{ average_x, average_y };

    std.mem.sort([2]f32, vertices, .{ .centroid = centroid }, Fns.cmpfn);
}

test "make_polygon" {
    const vertices: [][2]f32 = &.{ &.{ 0.0, 0.0 }, &.{ 1.0, 1.0 }, &.{ 1.0, 0.0 }, &.{ 0.0, 1.0 } };
    make_polygon(vertices);
    try std.testing.expect(vertices[0][0] == 0.0 and vertices[0][1] == 0.0);
    try std.testing.expect(vertices[1][0] == 1.0 and vertices[1][1] == 0.0);
    try std.testing.expect(vertices[2][0] == 1.0 and vertices[2][1] == 1.0);
    try std.testing.expect(vertices[3][0] == 0.0 and vertices[3][1] == 1.0);
}
