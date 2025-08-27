const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const clap_mod = b.addModule("clap", .{
        .root_source_file = b.path("src/clap.zig"),
        .target = target,
        .optimize = optimize,
    });
    _ = clap_mod;

    const draw_mod = b.addModule("draw", .{
        .root_source_file = b.path("src/draw.zig"),
        .target = target,
        .optimize = optimize,
    });
    _ = draw_mod;

    const io_adapter_mod = b.addModule("io_adapter", .{
        .root_source_file = b.path("src/io-adapter.zig"),
        .target = target,
        .optimize = optimize,
    });
    _ = io_adapter_mod;

    const misc_mod = b.addModule("misc", .{
        .root_source_file = b.path("src/misc.zig"),
        .target = target,
        .optimize = optimize,
    });
    _ = misc_mod;
}