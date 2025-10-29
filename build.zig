const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const clap_mod = b.addModule("clap", .{
        .root_source_file = b.path("src/clap.zig"),
        .target = target,
        .optimize = optimize,
    });

    const draw_mod = b.addModule("draw", .{
        .root_source_file = b.path("src/draw.zig"),
        .target = target,
        .optimize = optimize,
    });

    const io_adapter_mod = b.addModule("io_adapter", .{
        .root_source_file = b.path("src/io-adapter.zig"),
        .target = target,
        .optimize = optimize,
    });

    const misc_mod = b.addModule("misc", .{
        .root_source_file = b.path("src/misc.zig"),
        .target = target,
        .optimize = optimize,
    });

    const geometry_mod = b.addModule("geometry", .{
        .root_source_file = b.path("src/geometry.zig"),
        .target = target,
        .optimize = optimize,
    });

    const imports: []const std.Build.Module.Import = &.{
        .{ .name = "clap", .module = clap_mod },
        .{ .name = "draw", .module = draw_mod },
        .{ .name = "ioAdapter", .module = io_adapter_mod },
        .{ .name = "misc", .module = misc_mod },
        .{ .name = "geometry", .module = geometry_mod },
    };

    add_example(b, "draw", "examples/draw.zig", "Run the draw example", imports, target, optimize);
    add_test(b, "draw-test", draw_mod);
    add_test(b, "geometry-test", geometry_mod);
}

pub fn add_test(b: *std.Build, command: []const u8, mod: *std.Build.Module) void {
    const exe_unit_tests = b.addTest(.{
        .root_module = mod,
    });

    const run_exe_unit_tests = b.addRunArtifact(exe_unit_tests);

    // Similar to creating the run step earlier, this exposes a `test` step to
    // the `zig build --help` menu, providing a way for the user to request
    // running the unit tests.
    const test_step = b.step(command, "Run unit tests");
    test_step.dependOn(&run_exe_unit_tests.step);
}

pub fn add_example(b: *std.Build, command: []const u8, source: []const u8,
    help: []const u8, imports: []const std.Build.Module.Import,
    target: std.Build.ResolvedTarget, optimize: std.builtin.OptimizeMode) void {
    const exe_mod = b.addExecutable(.{
        .name = command,
        .root_module = b.createModule(.{
            .root_source_file = b.path(source),
            .target = target,
            .optimize = optimize,
            .imports = imports,
        }),
    });
    exe_mod.linkSystemLibrary("SDL2");
    exe_mod.linkSystemLibrary("c");
    b.installArtifact(exe_mod);
    const run_cmd = b.addRunArtifact(exe_mod);
    run_cmd.step.dependOn(b.getInstallStep());
    if (b.args) |args| {
        run_cmd.addArgs(args);
    }

    const run_step = b.step(command, help);
    run_step.dependOn(&run_cmd.step);
}
