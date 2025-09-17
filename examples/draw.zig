const std = @import("std");
const draw = @import("draw");
const ioAdapter = @import("ioAdapter");

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

    if (std.mem.eql(u8, example, "square")) {
        try runSquareExample();
    } else {
        std.debug.print("Unknown example: {s}\n", .{example});
    }
}

fn runSquareExample() !void {
    const width = 800;
    const height = 600;
    const num_squares = 100;

    var adapter = try ioAdapter.SDLAdapter.init(width, height);
    defer adapter.deinit();

    const ctx = draw.DrawContext(width, height);

    // Use a fixed seed for repeatable randomness
    var prng = std.rand.DefaultPrng.init(12345);

    for (num_squares) |_| {
        // Random position and size
        const x = prng.random().intRangeAtMost(u32, 0, width - 50);
        const y = prng.random().intRangeAtMost(u32, 0, height - 50);
        const size = prng.random().intRangeAtMost(u32, 20, 80);

        // Random color
        const r = prng.random().intRangeAtMost(u8, 0, 255);
        const g = prng.random().intRangeAtMost(u8, 0, 255);
        const b = prng.random().intRangeAtMost(u8, 0, 255);

        ctx.setStrokeColor(draw.Color{ .r = r, .g = g, .b = b, .a = 255 });
        ctx.strokeRect(x, y, size, size);
    }

    adapter.present();

    std.debug.print("Press any key to exit...\n", .{});
    _ = ioAdapter.waitForKey();
}
