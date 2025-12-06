// zig build -freference-trace=8 draw -- line
// zig build -freference-trace=8 draw -- fillPolygon
// zig build -freference-trace=8 draw -- png
const std = @import("std");
const misc = @import("misc");
const draw = @import("draw");
const ioAdapter = @import("ioAdapter");
const geometry = @import("geometry");
const png = @import("png");

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
    } else {
        std.debug.print("Unknown example: {s}\n", .{ example });
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
        _ = adapter.interface.waitForKey();
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
        _ = adapter.interface.waitForKey();
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
        _ = adapter.interface.waitForKey();
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