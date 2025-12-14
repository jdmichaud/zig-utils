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
    } else {
        std.debug.print("Unknown example: {s}\n", .{ example });
    }
}

fn runPrintExample(allocator: std.mem.Allocator, ttf_path: []const u8, text: []const u8) !void {
    const width = 800;
    const height = 600;

    var adapter = try ioAdapter.SDLAdapter.init(width, height);
    defer adapter.deinit();

    var context = try draw.DrawContext.init(allocator, width, height);
    context.alpha = true;
    context.imageSmoothingEnabled = true;
    context.imageSmoothingQuality = "low";

    // Load the PNG file
    const input = try misc.load(ttf_path);
    defer std.posix.munmap(input);
    var ttf_font = try ttf.TtfFont.load(allocator, input);
    defer ttf_font.deinit(allocator);

    var quit = false;
    var render_time: i128 = 0;

    const font_width = 100;
    const font_height = 100;

    var x: i32 = @intCast(@divTrunc(width, 2) - @divTrunc(font_width, 2));
    var y: i32 = @intCast(@divTrunc(height, 2) - @divTrunc(font_height, 2));
    const angle = 0;
    var zoom_factor: f32 = 1.0;

    var mouse_position: [2]i32 = .{ 0, 0 };
    var mouse_buttons: u3 = 0;
    _ = text;
    while (!quit) {
        misc.x11checkerboard(width, height, context.buffer);
        context.setTransform(1, 0, 0, 1, 0, 0);
        context.translate(misc.asf32(x) + misc.asf32(font_width) / 2, misc.asf32(y) + misc.asf32(font_height) / 2);
        context.rotate(angle);
        context.scale(zoom_factor, zoom_factor);

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

        // If a dot went out of the screen, just remove it and start the loop again
        vertices.clearRetainingCapacity();
        for (dots.items, 0..) |*dot, doti| {
            const newx = dot.x + dot.speedx;
            const newy = dot.y + dot.speedy;
            if (newx < 0 or newx >= width or newy < 0 or newy >= height) {
                _ = dots.swapRemove(doti);
                break;
            }
            dot.x = newx;
            dot.y = newy;
            try vertices.append(.{ newx, newy });
        }

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

