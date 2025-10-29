// zig build -freference-trace=8 draw -- line
const std = @import("std");
const draw = @import("draw");
const ioAdapter = @import("ioAdapter");
const geometry = @import("geometry");

const stdout = std.io.getStdOut().writer();
const stderr = std.io.getStdErr().writer();

const Shape = enum {
    line,
    rect,
    fillrect,
    fillcircle
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
    } else {
        std.debug.print("Unknown example: {s}\n", .{ example });
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