// zig build -freference-trace=8 draw -- line
// zig build -freference-trace=8 draw -- fillPolygon
// zig build -freference-trace=8 draw -- png examples/sprite.png
// zig build -freference-trace=8 draw -- draw examples/sprite.png
// zig build -freference-trace=8 draw -- interaction examples/sprite.png
// zig build -freference-trace=8 draw -- print /usr/share/fonts/truetype/dejavu/DejaVuSans.ttf Test
const std = @import("std");
const misc = @import("misc");
const draw = @import("draw");
const ioAdapter = @import("io_adapter");
const geometry = @import("geometry");
const png = @import("png");
const ttf = @import("ttf");
const zlm = @import("zlm");

const stdout = std.io.getStdOut().writer();
const stderr = std.io.getStdErr().writer();

const Shape = enum {
    line,
    rect,
    fillrect,
    fillcircle,
};

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    const allocator = gpa.allocator();

    const args = try std.process.argsAlloc(allocator);
    defer std.process.argsFree(allocator, args);

    if (args.len < 2) {
        std.debug.print("Usage: {s} <example>\n", .{args[0]});
        return;
    }

    const example = args[1];

    if (std.meta.stringToEnum(Shape, example)) |shape| {
        try runShapeExample(allocator, shape);
    } else if (std.mem.eql(u8, example, "strokePolygon")) {
        try runPolygonExample(allocator, .stroke);
    } else if (std.mem.eql(u8, example, "fillPolygon")) {
        try runPolygonExample(allocator, .fill);
    } else if (std.mem.eql(u8, example, "png")) {
        if (args.len != 3) {
            std.debug.print("Usage: {s} png <path/to/image.png>\n", .{args[0]});
            return;
        }
        try runPngExample(allocator, args[2]);
    } else if (std.mem.eql(u8, example, "draw")) {
        if (args.len != 3) {
            std.debug.print("Usage: {s} draw <path/to/image.png>\n", .{args[0]});
            return;
        }
        try runDrawExample(allocator, args[2]);
    } else if (std.mem.eql(u8, example, "interaction")) {
        if (args.len != 3) {
            std.debug.print("Usage: {s} interaction <path/to/image.png>\n", .{args[0]});
            return;
        }
        try runInteractionExample(allocator, args[2]);
    } else if (std.mem.eql(u8, example, "print")) {
        if (args.len != 4) {
            std.debug.print("Usage: {s} print <path/to/ttf_font.ttf> <text>\n", .{args[0]});
            return;
        }
        try runPrintExample(allocator, args[2], args[3]);
    } else if (std.mem.eql(u8, example, "paint")) {
        if (args.len != 2) {
            std.debug.print("Usage: {s} paint\n", .{args[0]});
            return;
        }
        try runPaintExample(allocator);
    } else {
        std.debug.print("Unknown example: {s}\n", .{ example });
    }
}

fn runPaintExample(allocator: std.mem.Allocator) !void {
    const screen_width = 800;
    const screen_height = 600;
    const canvas_width = screen_width;
    const canvas_height = screen_height;

    var adapter = try ioAdapter.SDLAdapter.init(screen_width, screen_height);
    defer adapter.deinit();

    var context = try draw.DrawContext.init(allocator, canvas_width, canvas_height);
    context.alpha = true;
    context.imageSmoothingEnabled = true;
    context.imageSmoothingQuality = "low";

    var quit = false;
    var render_time: i128 = 0;

    var mouse_position: [2]i16 = .{ 0, 0 };
    var mouse_buttons: u3 = 0;

    var x: f32 = 0;
    var y: f32 = 0;
    const angle = 0;
    var zoom_factor: f32 = 1.0;

    // State machine handling the creation of lines and curves
    // S to change from line to curves
    // Escape to stop the current line/curve
    const StateMachine = struct {
        const Line = struct { start_point: @Vector(2, i16), end_point: @Vector(2, i16) };
        const Curve = struct { start_point: @Vector(2, i16), control_point: @Vector(2, i16), end_point: @Vector(2, i16) };
        const LineState = enum { first_point, next_point };
        const CurveState = enum { first_point, control_point, second_point };
        const State = union(enum) {
            line: LineState,
            curve: CurveState,
        };

        state: State,
        lines: std.ArrayList(Line),
        curves: std.ArrayList(Curve),
        start_point: ?@Vector(2, i16) = null,
        control_point: ?@Vector(2, i16) = null,

        mouse_position: @Vector(2, i16) = .{ 0, 0 },

        const Self = @This();
        pub fn consume(self: *Self, event: ioAdapter.InputEvent) void {
            var mouse_buttons2: u3 = 0;
            switch (event) {
                ioAdapter.EventType.KeyDown => |keyEvent| {
                    switch (keyEvent.scancode) {
                        ioAdapter.Scancode.ESCAPE => {
                            self.start_point = null;
                            self.control_point = null;
                            switch (self.state) {
                                .line => self.state = .{ .line = .first_point },
                                .curve => self.state = .{ .curve = .first_point },
                            }
                        },
                        ioAdapter.Scancode.S => {
                            self.start_point = null;
                            self.control_point = null;
                            switch (self.state) {
                                .line => self.state = .{ .curve = .first_point },
                                .curve => self.state = .{ .line = .first_point },
                            }
                        },
                        else => {},
                    }
                },
                ioAdapter.EventType.MouseMove => |payload| {
                    self.mouse_position[0] = @intCast(payload.x);
                    self.mouse_position[1] = @intCast(payload.y);
                },
                ioAdapter.EventType.MouseDown => |payload| {
                    const which_btn = misc.asInt(u2, @intFromEnum(payload.button) - 1);
                    mouse_buttons2 |= @as(u3, 1) << which_btn;

                    if (mouse_buttons2 & @intFromEnum(ioAdapter.MouseButton.Left) != 0) {
                        switch (self.state) {
                            .line => |step| {
                                switch (step) {
                                    .first_point => {
                                        self.start_point = .{ self.mouse_position[0], self.mouse_position[1] };
                                        self.state = .{ .line = .next_point };
                                    },
                                    .next_point => {
                                        const current_point = .{ self.mouse_position[0], self.mouse_position[1] };
                                        self.lines.append(Line{
                                            .start_point = self.start_point.?,
                                            .end_point = current_point,
                                        }) catch {};
                                        self.start_point = current_point;
                                    },
                                }
                            },
                            .curve => |step| {
                                switch (step) {
                                    .first_point => {
                                        self.start_point = .{ self.mouse_position[0], self.mouse_position[1] };
                                        self.state = .{ .curve = .control_point };
                                    },
                                    .control_point => {
                                        self.control_point = .{ self.mouse_position[0], self.mouse_position[1] };
                                        self.state = .{ .curve = .second_point };
                                    },
                                    .second_point => {
                                        const current_point = .{ self.mouse_position[0], self.mouse_position[1] };
                                        self.curves.append(Curve{
                                            .start_point = self.start_point.?,
                                            .control_point = self.control_point.?,
                                            .end_point = current_point,
                                        }) catch {};
                                        self.start_point = null;
                                        self.control_point = null;
                                        self.state = .{ .curve = .first_point };
                                    },
                                }
                            },
                        }
                    }
                },
                else => {},
            }
        }
    };
    var state_machine: StateMachine = .{
        .state = .{ .line = .first_point },
        .lines = std.ArrayList(StateMachine.Line).init(allocator),
        .curves = std.ArrayList(StateMachine.Curve).init(allocator),
    };
    defer {
        state_machine.lines.deinit();
        state_machine.curves.deinit();
    }

    while (!quit) {
        misc.x11checkerboard(canvas_width, canvas_height, context.buffer);
        context.setTransform(1, 0, 0, 1, 0, 0);
        switch (state_machine.state) {
            .line => context.fillText("line", 10, 10),
            .curve => context.fillText("curve", 10, 10),
        }
        context.rotate(angle);
        context.translate(screen_width / 2, screen_height / 2);
        context.scale(zoom_factor, zoom_factor);
        context.translate(-screen_width / 2, -screen_height / 2);

        context.translate(x, y);

        for (state_machine.lines.items) |line| {
            context.strokeStyle = draw.DrawContext.shadeColor(255, 255, 255);
            context.moveTo(line.start_point[0], line.start_point[1]);
            context.lineTo(line.end_point[0], line.end_point[1]);
            context.stroke();
        }

        for (state_machine.curves.items) |curve| {
            context.strokeStyle = draw.DrawContext.shadeColor(160, 160, 160);
            context.moveTo(curve.start_point[0], curve.start_point[1]);
            context.lineTo(curve.control_point[0], curve.control_point[1]);
            context.moveTo(curve.control_point[0], curve.control_point[1]);
            context.lineTo(curve.end_point[0], curve.end_point[1]);
            context.stroke();

            context.strokeStyle = draw.DrawContext.shadeColor(255, 255, 255);
            context.moveTo(curve.start_point[0], curve.start_point[1]);
            context.quadraticCurveTo(curve.control_point[0], curve.control_point[1], curve.end_point[0], curve.end_point[1]);
            context.stroke();
        }

        if (state_machine.start_point) |start_point| {
            switch (state_machine.state) {
                .line => {
                    context.strokeStyle = draw.DrawContext.shadeColor(255, 255, 255);
                    context.moveTo(start_point[0], start_point[1]);
                    context.lineTo(state_machine.mouse_position[0], state_machine.mouse_position[1]);
                    context.stroke();
                },
                .curve => {
                    context.moveTo(start_point[0], start_point[1]);
                    if (state_machine.control_point) |control_point| {
                        context.strokeStyle = draw.DrawContext.shadeColor(160, 160, 160);
                        context.thickness = 4;
                        context.moveTo(start_point[0], start_point[1]);
                        context.lineTo(control_point[0], control_point[1]);
                        context.lineTo(mouse_position[0], mouse_position[1]);
                        context.stroke();
                        context.strokeStyle = draw.DrawContext.shadeColor(255, 255, 255);
                        context.thickness = 1;
                        context.moveTo(start_point[0], start_point[1]);
                        context.quadraticCurveTo(
                            control_point[0], control_point[1],
                            mouse_position[0], mouse_position[1],
                        );
                        context.stroke();
                    } else {
                        context.strokeStyle = draw.DrawContext.shadeColor(160, 160, 160);
                        context.thickness = 4;
                        context.moveTo(start_point[0], start_point[1]);
                        context.lineTo(mouse_position[0], mouse_position[1]);
                        context.stroke();
                        context.strokeStyle = draw.DrawContext.shadeColor(255, 160, 160);
                        context.thickness = 1;
                        context.fillCircle(mouse_position[0], mouse_position[1], 2);
                        context.stroke();
                    }
                },
            }
        }

        while (adapter.interface.getEvent()) |event| {
            switch (event) {
                ioAdapter.EventType.KeyDown => |keyEvent| {
                    switch (keyEvent.scancode) {
                        ioAdapter.Scancode.ESCAPE => switch (state_machine.state) {
                            .line => |step| switch (step) {
                                .first_point => quit = true,
                                else => {},
                            },
                            .curve => |step| switch (step) {
                                .first_point => quit = true,
                                .control_point => {},
                                .second_point => {},
                            },
                        },
                        else => {},
                    }
                },
                ioAdapter.EventType.MouseWheel => |payload| {
                    zoom_factor = @max(0.1, zoom_factor + if (payload.y < 0) @as(f32, -0.1) else @as(f32, 0.1));
                },
                ioAdapter.EventType.MouseDown => |payload| {
                    const which_btn = misc.asInt(u2, @intFromEnum(payload.button) - 1);
                    mouse_buttons |= @as(u3, 1) << which_btn;
                },
                ioAdapter.EventType.MouseUp => |payload| {
                    const which_btn = misc.asInt(u2, (@intFromEnum(payload.button) - 1));
                    mouse_buttons ^= mouse_buttons & (@as(u3, 1) << which_btn);
                },
                ioAdapter.EventType.MouseMove => |payload| {
                    mouse_position[0] = @intCast(payload.x);
                    mouse_position[1] = @intCast(payload.y);
                    if (mouse_buttons & (@as(u8, 1) << @intFromEnum(ioAdapter.MouseButton.Right) - 1) != 0) {
                        x += misc.asf32(payload.dx) / zoom_factor;
                        y += misc.asf32(payload.dy) / zoom_factor;
                    }
                },
                else => {},
            }
            state_machine.consume(event);
            // std.log.debug("state_machine {any}", .{ state_machine.state });
        }

        // Only go for 60fps
        const now = std.time.nanoTimestamp();
        if (quit or now - render_time > 16000000) {
            // context.clearRect(0, 0, context.width, context.height);
            adapter.interface.drawImage(context.buffer, 0, 0, canvas_width, canvas_height);
            adapter.interface.renderScene();
            render_time = now;
        }
    }
}

fn runPrintExample(allocator: std.mem.Allocator, ttf_path: []const u8, text: []const u8) !void {
    const screen_width = 800;
    const screen_height = 600;
    const canvas_width = screen_width;
    const canvas_height = screen_height;

    var adapter = try ioAdapter.SDLAdapter.init(screen_width, screen_height);
    defer adapter.deinit();

    var context = try draw.DrawContext.init(allocator, canvas_width, canvas_height);
    context.alpha = true;
    context.imageSmoothingEnabled = true;
    context.imageSmoothingQuality = "low";

    // Load the TTF file
    const input = try misc.load(ttf_path);
    defer std.posix.munmap(input);
    var ttf_font = try ttf.TtfFont.load(allocator, input);
    defer ttf_font.deinit(allocator);

    for (ttf_font.table_records) |tr| {
        std.log.debug("{s} 0x{x:0>8} {}", .{ tr.tag, tr.offset, tr.length });
    }

    ttf_font.maxp_table.pretty_print();
    ttf_font.head_table.pretty_print();
    ttf_font.cmap_table.pretty_print();
    std.log.debug("{} glyphs", .{ ttf_font.loca_table.count() });

    for (text) |c| {
        var glyph = try ttf_font.getGlyph(allocator, c);
        defer glyph.deinit(allocator);
        var contours = glyph.contours();
        std.log.debug("{c} {} contours", .{ c, glyph.number_of_contours });
        std.log.debug("  bbox {d} {d} {d} {d}", .{ glyph.x_min, glyph.y_min, glyph.x_max, glyph.y_max });
        var i: usize = 0;
        while (contours.next()) |points| {
            std.log.debug("contour {}", .{ i });
            for (points) |point| {
                std.log.debug("  {any}", .{ point });
            }
            i += 1;
        }
    }

    var quit = false;
    var render_time: i128 = 0;

    var x: f32 = 0;
    var y: f32 = 0;
    const angle = 0;
    var zoom_factor: f32 = 1.0;

    var mouse_position: [2]i32 = .{ 0, 0 };
    var mouse_buttons: u3 = 0;

    while (!quit) {
        misc.x11checkerboard(canvas_width, canvas_height, context.buffer);
        context.setTransform(1, 0, 0, 1, 0, 0);
        context.rotate(angle);
        context.translate(screen_width / 2, screen_height / 2);
        context.scale(zoom_factor, zoom_factor);
        context.translate(-screen_width / 2, -screen_height / 2);

        context.translate(x, y);

        context.strokeStyle = 0xFFBBBBBB;
        context.fillStyle = 0xFF0000FF;

        var penx: i32 = 0;
        for (text) |c| {
            var glyph = try ttf_font.getGlyph(allocator, c);
            const metric = try ttf_font.getGlyphMetric(c);
            defer glyph.deinit(allocator);
            const t_height = ttf_font.head_table.y_max - ttf_font.head_table.y_min;
            const t_ratio = misc.asf32(canvas_height) / misc.asf32(t_height) * 0.5;

            const left_bound = penx; // + metric.lsb;
            var contours = glyph.contours();
            while (contours.next()) |points| {
                var first = true;
                context.beginPath();
                for (points) |point| {
                    const point_x = misc.asInt(i16, misc.asf32(point.x + left_bound) * t_ratio);
                    const point_y = misc.asInt(i16, canvas_height - misc.asf32(point.y) * t_ratio);
                    if (point.on_curve) {
                        if (first) {
                            first = false;
                            context.moveTo(point_x, point_y);
                        } else {
                            context.lineTo(point_x, point_y);
                        }
                    } else {
                        context.fillCircle(point_x, point_y, 3);
                    }
                }
                context.closePath();
                context.stroke();
            }

            var curves = glyph.curves();
            context.strokeStyle = 0xFFFF0000;
            context.beginPath();
            while (curves.next()) |segment| {
                switch (segment) {
                    .move_to => |point| {
                        const point_x = misc.asInt(i16, misc.asf32(point.x + left_bound) * t_ratio);
                        const point_y = misc.asInt(i16, canvas_height - misc.asf32(point.y) * t_ratio);
                        context.moveTo(point_x, point_y);
                    },
                    .line_to => |point| {
                        const point_x = misc.asInt(i16, misc.asf32(point.x + left_bound) * t_ratio);
                        const point_y = misc.asInt(i16, canvas_height - misc.asf32(point.y) * t_ratio);
                        context.lineTo(point_x, point_y);
                    },
                    .quad_to => |parameters| {
                        const point_x = misc.asInt(i16, misc.asf32(parameters.x + left_bound) * t_ratio);
                        const point_y = misc.asInt(i16, canvas_height - misc.asf32(parameters.y) * t_ratio);
                        const control_point_x = misc.asInt(i16, misc.asf32(parameters.cx + left_bound) * t_ratio);
                        const control_point_y = misc.asInt(i16, canvas_height - misc.asf32(parameters.cy) * t_ratio);
                        context.quadraticCurveTo(control_point_x, control_point_y, point_x, point_y);
                    },
                }
            }
            context.closePath();
            context.stroke();
            penx += metric.advance_width;
        }

        while (adapter.interface.getEvent()) |event| {
            switch (event) {
                ioAdapter.EventType.KeyDown => |keyEvent| {
                    switch (keyEvent.scancode) {
                        ioAdapter.Scancode.ESCAPE => quit = true,
                        ioAdapter.Scancode.I => {
                            if (std.mem.eql(u8, context.imageSmoothingQuality, "low")) {
                                context.imageSmoothingQuality = "medium";
                            } else if (std.mem.eql(u8, context.imageSmoothingQuality, "medium")) {
                                context.imageSmoothingQuality = "high";
                            } else if (std.mem.eql(u8, context.imageSmoothingQuality, "high")) {
                                context.imageSmoothingQuality = "low";
                            }
                            try stdout.print("set interpolation to {s}\n", .{ context.imageSmoothingQuality });
                        },
                        else => {},
                    }
                },
                ioAdapter.EventType.MouseWheel => |payload| {
                    zoom_factor = @max(0.1, zoom_factor + if (payload.y < 0) @as(f32, -0.1) else @as(f32, 0.1));
                },
                ioAdapter.EventType.MouseDown => |payload| {
                    const which_btn = misc.asInt(u2, @intFromEnum(payload.button) - 1);
                    mouse_buttons |= @as(u3, 1) << which_btn;
                },
                ioAdapter.EventType.MouseUp => |payload| {
                    const which_btn = misc.asInt(u2, (@intFromEnum(payload.button) - 1));
                    mouse_buttons ^= mouse_buttons & (@as(u3, 1) << which_btn);
                },
                ioAdapter.EventType.MouseMove => |payload| {
                    mouse_position[0] = payload.x;
                    mouse_position[1] = payload.y;
                    if (mouse_buttons & @intFromEnum(ioAdapter.MouseButton.Left) != 0) {
                        x += misc.asf32(payload.dx) / zoom_factor;
                        y += misc.asf32(payload.dy) / zoom_factor;
                    }
                },
                else => {},
            }
        }

        // Only go for 60fps
        const now = std.time.nanoTimestamp();
        if (quit or now - render_time > 16000000) {
            // context.clearRect(0, 0, context.width, context.height);
            adapter.interface.drawImage(context.buffer, 0, 0, canvas_width, canvas_height);
            adapter.interface.renderScene();
            render_time = now;
        }
    }
}

fn runInteractionExample(allocator: std.mem.Allocator, image_path: []const u8) !void {
    const width = 800;
    const height = 600;

    var adapter = try ioAdapter.SDLAdapter.init(width, height);
    defer adapter.deinit();

    var context = try draw.DrawContext.init(allocator, width, height);
    context.alpha = true;
    context.imageSmoothingEnabled = true;
    context.imageSmoothingQuality = "low";

    // Load the PNG file
    const input = try misc.load(image_path);
    defer std.posix.munmap(input);
    var png_image_data = try png.load_png(allocator, input);
    defer png_image_data.deinit(allocator);

    const imageData = draw.ImageData{
        .width = png_image_data.header.getWidth(),
        .height = png_image_data.header.getHeight(),
        .pixel_format = draw.PixelFormat.RGBA8,
        .data = png_image_data.data,
    };

    var quit = false;
    var render_time: i128 = 0;

    var x: i32 = @intCast(@divTrunc(width, 2) - @divTrunc(imageData.width, 2));
    var y: i32 = @intCast(@divTrunc(height, 2) - @divTrunc(imageData.height, 2));
    const angle = 0;
    var zoom_factor: f32 = 1.0;

    var mouse_position: [2]i32 = .{ 0, 0 };
    var mouse_buttons: u3 = 0;
    var shift_pressed = false;
    var ctrl_pressed = false;
    while (!quit) {
        misc.x11checkerboard(width, height, context.buffer);
        context.setTransform(1, 0, 0, 1, 0, 0);
        context.translate(misc.asf32(x) + misc.asf32(imageData.width) / 2, misc.asf32(y) + misc.asf32(imageData.height) / 2);
        context.rotate(angle);
        context.scale(zoom_factor, zoom_factor);
        context.drawImage(imageData, -@as(i32, @intCast(imageData.width / 2)), -@as(i32, @intCast(imageData.height / 2)));

        while (adapter.interface.getEvent()) |event| {
            switch (event) {
                ioAdapter.EventType.KeyDown => |keyEvent| {
                    switch (keyEvent.scancode) {
                        ioAdapter.Scancode.ESCAPE => quit = true,
                        ioAdapter.Scancode.LSHIFT => shift_pressed = true,
                        ioAdapter.Scancode.RSHIFT => shift_pressed = true,
                        ioAdapter.Scancode.RCTRL => ctrl_pressed = true,
                        ioAdapter.Scancode.LCTRL => ctrl_pressed = true,
                        ioAdapter.Scancode.I => {
                            if (std.mem.eql(u8, context.imageSmoothingQuality, "low")) {
                                context.imageSmoothingQuality = "medium";
                            } else if (std.mem.eql(u8, context.imageSmoothingQuality, "medium")) {
                                context.imageSmoothingQuality = "high";
                            } else if (std.mem.eql(u8, context.imageSmoothingQuality, "high")) {
                                context.imageSmoothingQuality = "low";
                            }
                            try stdout.print("set interpolation to {s}\n", .{ context.imageSmoothingQuality });
                        },
                        else => {},
                    }
                },
                ioAdapter.EventType.KeyUp => |keyEvent| {
                    switch (keyEvent.scancode) {
                        ioAdapter.Scancode.LSHIFT => shift_pressed = false,
                        ioAdapter.Scancode.RSHIFT => shift_pressed = false,
                        ioAdapter.Scancode.RCTRL => ctrl_pressed = false,
                        ioAdapter.Scancode.LCTRL => ctrl_pressed = false,
                        else => {},
                    }
                },
                ioAdapter.EventType.MouseWheel => |payload| {
                    const factor: f32 = if (shift_pressed) 0.5 else 0.1;
                    zoom_factor = @max(0.1, zoom_factor + if (payload.y < 0) @as(f32, -factor) else @as(f32, factor));
                },
                ioAdapter.EventType.MouseDown => |payload| {
                    const which_btn = misc.asInt(u2, @intFromEnum(payload.button) - 1);
                    mouse_buttons |= @as(u3, 1) << which_btn;
                },
                ioAdapter.EventType.MouseUp => |payload| {
                    const which_btn = misc.asInt(u2, (@intFromEnum(payload.button) - 1));
                    mouse_buttons ^= mouse_buttons & (@as(u3, 1) << which_btn);
                },
                ioAdapter.EventType.MouseMove => |payload| {
                    mouse_position[0] = payload.x;
                    mouse_position[1] = payload.y;
                    if (mouse_buttons & @intFromEnum(ioAdapter.MouseButton.Left) != 0) {
                        x += payload.dx;
                        y += payload.dy;
                    }
                },
                else => {},
            }
        }

        // Only go for 60fps
        const now = std.time.nanoTimestamp();
        if (quit or now - render_time > 16000000) {
            // context.clearRect(0, 0, context.width, context.height);
            adapter.interface.drawImage(context.buffer, 0, 0, width, height);
            adapter.interface.renderScene();
            render_time = now;
        }
    }
}

fn runDrawExample(allocator: std.mem.Allocator, image_path: []const u8) !void {
    const width = 800;
    const height = 600;
    const num_image = 100000;
    // const num_image = 1;

    var adapter = try ioAdapter.SDLAdapter.init(width, height);
    defer adapter.deinit();

    var context = try draw.DrawContext.init(allocator, width, height);
    context.alpha = true;
    // context.imageSmoothingEnabled = true;
    // context.imageSmoothingQuality = "high";
    // to check the transparency
    misc.x11checkerboard(width, height, context.buffer);

    // Use a fixed seed for repeatable randomness
    var prng = std.Random.DefaultPrng.init(12345);

    // Load the PNG file
    const input = try misc.load(image_path);
    defer std.posix.munmap(input);
    var png_image_data = try png.load_png(allocator, input);
    defer png_image_data.deinit(allocator);

    const imageData = draw.ImageData{
        .width = png_image_data.header.getWidth(),
        .height = png_image_data.header.getHeight(),
        .pixel_format = draw.PixelFormat.RGBA8,
        .data = png_image_data.data,
    };

    var quit = false;
    var i: u32 = 0;
    var render_time: i128 = 0;
    const start_time = std.time.nanoTimestamp();
    var ticks: usize = 0;
    while (i < num_image and !quit) {
        i += 1;
        // Random position and size
        const x = prng.random().intRangeAtMost(i16, -@as(i16, @intCast(imageData.width)), width);
        const y = prng.random().intRangeAtMost(i16, -@as(i16, @intCast(imageData.height)), height);
        const angle = prng.random().float(f32) * std.math.pi;
        const zoom_factor = prng.random().float(f32) * 1.4 + 0.01;
        // _ = prng;
        // const x = 100;
        // const y = 100;
        // const angle = std.math.pi / 2.0;
        // const zoom_factor = 2;

        context.setTransform(1, 0, 0, 1, 0, 0);
        context.translate(misc.asf32(x) + misc.asf32(imageData.width) / 2, misc.asf32(y) + misc.asf32(imageData.height) / 2);
        context.rotate(angle);
        context.scale(zoom_factor, zoom_factor);
        context.drawImage(imageData, -@as(i32, @intCast(imageData.width / 2)), -@as(i32, @intCast(imageData.height / 2)));

        quit = isQuit(&adapter.interface);

        ticks += 1;
        // Only go for 60fps
        const now = std.time.nanoTimestamp();
        if (i >= num_image or quit or now - render_time > 16000000) {
            adapter.interface.drawImage(context.buffer, 0, 0, width, height);
            adapter.interface.renderScene();
            render_time = now;
            try stdout.print("{} blits/frames\n", .{ ticks });
            ticks = 0;
        }
    }

    try stdout.print("Done in {}ns...\n", .{ @divTrunc(std.time.nanoTimestamp() - start_time, 1000000) });

    if (!quit) {
        try stdout.print("Press any key to exit...\n", .{});
        _ = waitForKey(&adapter.interface);
    }
}

fn runPngExample(allocator: std.mem.Allocator, image_path: []const u8) !void {
    const width = 800;
    const height = 600;
    const num_image = 100000;

    var adapter = try ioAdapter.SDLAdapter.init(width, height);
    defer adapter.deinit();

    var context = try draw.DrawContext.init(allocator, width, height);
    context.alpha = true;
    // to check the transparency
    misc.x11checkerboard(width, height, context.buffer);

    // Use a fixed seed for repeatable randomness
    var prng = std.Random.DefaultPrng.init(12345);

    // Load the PNG file
    const input = try misc.load(image_path);
    defer std.posix.munmap(input);
    var png_image_data = try png.load_png(allocator, input);
    defer png_image_data.deinit(allocator);

    const imageData = draw.ImageData{
        .width = png_image_data.header.getWidth(),
        .height = png_image_data.header.getHeight(),
        .pixel_format = draw.PixelFormat.RGBA8,
        .data = png_image_data.data,
    };

    var quit = false;
    var i: u32 = 0;
    var render_time: i128 = 0;
    var ticks: usize = 0;
    while (i < num_image and !quit) {
        i += 1;
        // Random position and size
        // const x = prng.random().intRangeAtMost(i16, 0, width - @as(i16, @intCast(imageData.width)));
        // const y = prng.random().intRangeAtMost(i16, 0, height - @as(i16, @intCast(imageData.height)));
        const x = prng.random().intRangeAtMost(i16, -@as(i16, @intCast(imageData.width)), width);
        const y = prng.random().intRangeAtMost(i16, -@as(i16, @intCast(imageData.height)), height);

        context.putImageData(imageData, x, y);

        quit = isQuit(&adapter.interface);

        ticks += 1;
        // Only go for 60fps
        const now = std.time.nanoTimestamp();
        if (i >= num_image or quit or now - render_time > 16000000) {
            adapter.interface.drawImage(context.buffer, 0, 0, width, height);
            adapter.interface.renderScene();
            render_time = now;
            try stdout.print("{} blits/frames\n", .{ ticks });
            ticks = 0;
        }
    }

    if (!quit) {
        try stdout.print("Press any key to exit...\n", .{});
        _ = waitForKey(&adapter.interface);
    }
}

fn runPolygonExample(allocator: std.mem.Allocator, action_type: enum { stroke, fill } ) !void {
    const width = 800;
    const height = 600;
    const num_dots = 4;

    var adapter = try ioAdapter.SDLAdapter.init(width, height);
    defer adapter.deinit();

    var context = try draw.DrawContext.init(allocator, width, height);

    // Use a fixed seed for repeatable randomness
    var prng = std.Random.DefaultPrng.init(12345);

    var quit = false;
    var ticks: usize = 0;
    const Dot = struct {
        x: f32,
        y: f32,
        speedx: f32,
        speedy: f32,
    };
    var dots = std.ArrayList(Dot).init(allocator);
    var vertices = std.ArrayList([2]f32).init(allocator);
    defer dots.deinit();
    while (!quit) {
        const then = std.time.nanoTimestamp();
        while (dots.items.len < num_dots) {
            const x = prng.random().float(f32) * (width - 300) + 100;
            const y = prng.random().float(f32) * (height - 300) + 100;
            const speedx = prng.random().float(f32) - 0.5;
            const speedy = prng.random().float(f32) - 0.5;
            try dots.append(.{ .x = x, .y = y, .speedx = speedx, .speedy = speedy });
        }

        // // If a dot went out of the screen, just remove it and start the loop again
        // vertices.clearRetainingCapacity();
        // std.log.debug("----", .{ });
        // for (dots.items, 0..) |*dot, doti| {
        //     const newx = dot.x + dot.speedx;
        //     const newy = dot.y + dot.speedy;
        //     if (newx < 0 or newx >= width or newy < 0 or newy >= height) {
        //         _ = dots.swapRemove(doti);
        //         break;
        //     }
        //     dot.x = newx;
        //     dot.y = newy;
        //     std.log.debug("{d} {d}", .{ newx, newy });
        //     try vertices.append(.{ newx, newy });
        // }

        // try vertices.append(.{ 453.65323, 45.756447 });
        // try vertices.append(.{ 135.27386, 366.7728 });
        // try vertices.append(.{ 282.59277, 259.66742 });
        // try vertices.append(.{ 295.758, 217.83728 });

        try vertices.append(.{ 4, 0 });
        try vertices.append(.{ 1, 3 });
        try vertices.append(.{ 3, 3 });
        try vertices.append(.{ 3, 2 });


        geometry.make_polygon(vertices.items);

        @memset(context.buffer, 0);
        context.strokeStyle = 0xFFBBBBBB;
        context.beginPath();
        for (vertices.items, 0..) |*vertex, vertexi| {
            const posx: i16 = @intFromFloat(@round(vertex[0]));
            const posy: i16 = @intFromFloat(@round(vertex[1]));
            context.fillCircle(posx, posy, 3);

            if (vertexi == 0) {
                context.moveTo(posx, posy);
            } else {
                context.lineTo(posx, posy);
            }
        }
        context.closePath();
        switch (action_type) {
            .stroke => context.stroke(),
            .fill => context.fill(),
        }

        adapter.interface.drawImage(context.buffer, 0, 0, width, height);
        adapter.interface.renderScene();

        ticks += 1;

        // Only go for 60fps
        const timePerFrame = std.time.nanoTimestamp() - then;
        if (timePerFrame < 16000000) {
            std.time.sleep(@intCast(16000000 - timePerFrame));
        }

        if (ticks % 100 == 0) {
            try stdout.print("{}ms\n", .{ @as(u16, @intFromFloat(@as(f32, @floatFromInt(timePerFrame)) / 1000000.0)) });
        }

        quit = isQuit(&adapter.interface);
    }

    if (!quit) {
        try stdout.print("Press any key to exit...\n", .{});
        _ = waitForKey(&adapter.interface);
    }
}

fn runShapeExample(allocator: std.mem.Allocator, shape: Shape) !void {
    const width = 800;
    const height = 600;
    const num_lines = 10000;

    var adapter = try ioAdapter.SDLAdapter.init(width, height);
    defer adapter.deinit();

    var context = try draw.DrawContext.init(allocator, width, height);

    // Use a fixed seed for repeatable randomness
    var prng = std.Random.DefaultPrng.init(12345);

    var quit = false;
    var i: u32 = 0;
    while (i < num_lines and !quit) {
        i += 1;
        // Random position and size
        const x = prng.random().intRangeAtMost(i16, 0, width);
        const y = prng.random().intRangeAtMost(i16, 0, height);
        const dx = prng.random().intRangeAtMost(i16, 0, width);
        const dy = prng.random().intRangeAtMost(i16, 0, height);

        // Random color
        const r = prng.random().intRangeAtMost(u8, 0, 255);
        const g = prng.random().intRangeAtMost(u8, 0, 255);
        const b = prng.random().intRangeAtMost(u8, 0, 255);

        context.fillStyle = draw.DrawContext.shadeColor(r, g, b);
        context.strokeStyle = draw.DrawContext.shadeColor(r, g, b);

        switch (shape) {
            .line => context.line(x, y, dx, dy),
            .rect => context.strokeRect(x, y, dx, dy),
            .fillrect => context.fillRect(x, y, dx, dy),
            .fillcircle => {
                const radius = prng.random().intRangeAtMost(i16, 0, 100);
                context.fillCircle(x, y, radius);
            }
        }

        var buffer: [256]u8 = undefined;
        const str = try std.fmt.bufPrintZ(&buffer, "{}", .{ i });
        context.fillText(str, 10, 10);
        adapter.interface.drawImage(context.buffer, 0, 0, width, height);
        adapter.interface.renderScene();

        quit = isQuit(&adapter.interface);
    }

    if (!quit) {
        try stdout.print("Press any key to exit...\n", .{});
        _ = waitForKey(&adapter.interface);
    }
}


fn isQuit(adapter: *ioAdapter.IOAdapter) bool {
    while (adapter.getEvent()) |event| {
        switch (event) {
            ioAdapter.EventType.KeyDown => |keyEvent| {
              if (keyEvent.scancode == ioAdapter.Scancode.ESCAPE)
                return true;
            },
            else => {},
        }
    }
    return false;
}

fn waitForKey(adapter: *ioAdapter.IOAdapter) void {
    loop: switch (adapter.waitEvent()) {
        .KeyDown => return,
        else => continue :loop adapter.waitEvent(),
    }
}

