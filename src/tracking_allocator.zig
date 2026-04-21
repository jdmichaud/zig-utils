//! Tracking allocator with per-call-site memory statistics.
//!
//! Wraps any other allocator and records, for each unique call site that
//! triggers alloc/resize/free: total bytes allocated, total bytes freed,
//! currently-live bytes, peak live bytes, and allocation count.
//!
//! Two capture modes, selected at init via `Config`:
//!
//!   • **Single-frame (default)**: records only the immediate return address
//!     that the allocator interface passes in. Cheap — one hash insert per
//!     alloc. Two paths that go through a common helper (like
//!     `std.ArrayList.append`) will be aggregated into the helper's site.
//!
//!   • **Full trace**: walks up to `MAX_FRAMES` frames using
//!     `std.debug.captureStackTrace`. Distinguishes callers even when they
//!     route through a common helper. More overhead — a stack walk per alloc.
//!
//! Usage:
//!     var gpa = std.heap.GeneralPurposeAllocator(.{}){};   // bookkeeping
//!     var tracker = TrackingAllocator.init(
//!         std.heap.page_allocator,                          // parent
//!         gpa.allocator(),                                  // backing
//!         .{ .capture_full_trace = false },                 // config
//!     );
//!     defer tracker.deinit();
//!
//!     const alloc = tracker.allocator();
//!     // ... use `alloc` everywhere in the program ...
//!
//!     try tracker.report(std.io.getStdErr().writer());
//!
//! Caveats:
//!   - Symbol resolution requires DWARF debug info (Debug/ReleaseSafe auto-ok).
//!   - Full-trace mode requires frame pointers or DWARF for the stack walk.
//!   - Single-threaded only; add a mutex for multi-threaded profiling.
//!   - `backing` MUST be a different allocator than `parent`, otherwise the
//!     bookkeeping's own allocations would show up in the profile and
//!     possibly cause recursion.
//!
//! Made with Claude Opus 4.7

const std = @import("std");

pub const Config = struct {
    /// When false (default), only the immediate return address is captured.
    /// When true, up to `MAX_FRAMES` call-stack frames are captured per
    /// allocation. The stack walk costs ~100ns-1µs per alloc.
    capture_full_trace: bool = false,
};

pub const TrackingAllocator = struct {
    parent: std.mem.Allocator,
    backing: std.mem.Allocator,
    config: Config,
    /// site_key → aggregated stats. Keyed by hash of the call trace.
    per_site: std.AutoHashMap(u64, SiteStats),
    /// live ptr → (size, site_key), so free can subtract correctly.
    live: std.AutoHashMap(usize, LiveAlloc),

    pub const MAX_FRAMES = 8;

    const LiveAlloc = struct {
        size: usize,
        site_key: u64,
    };

    const SiteStats = struct {
        count: usize = 0,
        bytes_allocated: usize = 0,
        bytes_freed: usize = 0,
        bytes_live: usize = 0,
        peak_live: usize = 0,
        /// Always populated; `trace_len` is 1 in single-frame mode.
        trace: [MAX_FRAMES]usize = undefined,
        trace_len: u8 = 0,
    };

    pub fn init(
        parent: std.mem.Allocator,
        backing: std.mem.Allocator,
        config: Config,
    ) TrackingAllocator {
        return .{
            .parent = parent,
            .backing = backing,
            .config = config,
            .per_site = std.AutoHashMap(u64, SiteStats).init(backing),
            .live = std.AutoHashMap(usize, LiveAlloc).init(backing),
        };
    }

    pub fn deinit(self: *TrackingAllocator) void {
        self.per_site.deinit();
        self.live.deinit();
    }

    pub fn allocator(self: *TrackingAllocator) std.mem.Allocator {
        return .{ .ptr = self, .vtable = &vtable };
    }

    const vtable: std.mem.Allocator.VTable = .{
        .alloc = alloc,
        .resize = resize,
        .remap = remap,
        .free = free,
    };

    fn captureTrace(self: *const TrackingAllocator, ret: usize) struct {
        addrs: [MAX_FRAMES]usize,
        len: u8,
    } {
        var addrs: [MAX_FRAMES]usize = @splat(0);
        var len: u8 = 0;

        if (self.config.capture_full_trace) {
            var st = std.builtin.StackTrace{
                .instruction_addresses = &addrs,
                .index = 0,
            };
            std.debug.captureStackTrace(ret, &st);
            len = @intCast(st.index);
        } else {
            addrs[0] = ret;
            len = 1;
        }

        return .{ .addrs = addrs, .len = len };
    }

    fn alloc(ctx: *anyopaque, n: usize, a: std.mem.Alignment, ret: usize) ?[*]u8 {
        const self: *TrackingAllocator = @ptrCast(@alignCast(ctx));
        const p = self.parent.rawAlloc(n, a, ret) orelse return null;

        const t = self.captureTrace(ret);
        const site_key = hashTrace(t.addrs[0..t.len]);
        self.recordAlloc(p, n, site_key, t.addrs, t.len);
        return p;
    }

    fn resize(ctx: *anyopaque, buf: []u8, a: std.mem.Alignment, n: usize, ret: usize) bool {
        const self: *TrackingAllocator = @ptrCast(@alignCast(ctx));
        if (!self.parent.rawResize(buf, a, n, ret)) return false;

        if (self.live.getPtr(@intFromPtr(buf.ptr))) |la| {
            const delta: isize = @as(isize, @intCast(n)) - @as(isize, @intCast(buf.len));
            la.size = n;
            if (self.per_site.getPtr(la.site_key)) |s| {
                self.applyDelta(s, delta);
            }
        }
        return true;
    }

    fn remap(ctx: *anyopaque, buf: []u8, a: std.mem.Alignment, n: usize, ret: usize) ?[*]u8 {
        const self: *TrackingAllocator = @ptrCast(@alignCast(ctx));
        const new_p = self.parent.rawRemap(buf, a, n, ret) orelse return null;

        // Move the live record from old ptr to new ptr.
        const old_entry = self.live.fetchRemove(@intFromPtr(buf.ptr));
        const site_key = if (old_entry) |oe| oe.value.site_key else blk: {
            // No prior record (shouldn't normally happen); capture fresh.
            const t = self.captureTrace(ret);
            const k = hashTrace(t.addrs[0..t.len]);
            const e = self.per_site.getOrPut(k) catch break :blk k;
            if (!e.found_existing) {
                e.value_ptr.* = .{ .trace = t.addrs, .trace_len = t.len };
            }
            break :blk k;
        };

        self.live.put(@intFromPtr(new_p), .{ .size = n, .site_key = site_key }) catch {};

        if (self.per_site.getPtr(site_key)) |s| {
            const delta: isize = @as(isize, @intCast(n)) - @as(isize, @intCast(buf.len));
            self.applyDelta(s, delta);
        }
        return new_p;
    }

    fn free(ctx: *anyopaque, buf: []u8, a: std.mem.Alignment, ret: usize) void {
        const self: *TrackingAllocator = @ptrCast(@alignCast(ctx));
        self.parent.rawFree(buf, a, ret);

        if (self.live.fetchRemove(@intFromPtr(buf.ptr))) |kv| {
            if (self.per_site.getPtr(kv.value.site_key)) |s| {
                s.bytes_freed += kv.value.size;
                s.bytes_live -= kv.value.size;
            }
        }
    }

    fn recordAlloc(
        self: *TrackingAllocator,
        p: [*]u8,
        n: usize,
        site_key: u64,
        trace: [MAX_FRAMES]usize,
        trace_len: u8,
    ) void {
        // Tolerate OOM in bookkeeping — the real alloc should still succeed.
        self.live.put(@intFromPtr(p), .{ .size = n, .site_key = site_key }) catch return;

        const e = self.per_site.getOrPut(site_key) catch return;
        if (!e.found_existing) {
            e.value_ptr.* = .{ .trace = trace, .trace_len = trace_len };
        }
        e.value_ptr.count += 1;
        e.value_ptr.bytes_allocated += n;
        e.value_ptr.bytes_live += n;
        if (e.value_ptr.bytes_live > e.value_ptr.peak_live) {
            e.value_ptr.peak_live = e.value_ptr.bytes_live;
        }
    }

    fn applyDelta(_: *TrackingAllocator, s: *SiteStats, delta: isize) void {
        if (delta > 0) {
            s.bytes_allocated += @intCast(delta);
            s.bytes_live += @intCast(delta);
        } else if (delta < 0) {
            s.bytes_freed += @intCast(-delta);
            s.bytes_live -= @intCast(-delta);
        }
        if (s.bytes_live > s.peak_live) s.peak_live = s.bytes_live;
    }

    fn hashTrace(trace: []const usize) u64 {
        return std.hash.Wyhash.hash(0, std.mem.sliceAsBytes(trace));
    }

    /// Total bytes currently held by all tracked allocations.
    pub fn totalLive(self: *const TrackingAllocator) usize {
        var sum: usize = 0;
        var it = self.per_site.valueIterator();
        while (it.next()) |v| sum += v.bytes_live;
        return sum;
    }

    /// Write a human-readable report to `writer`, sorted by peak live bytes
    /// descending. Each call site prints its single location (single-frame
    /// mode) or full stack trace (full-trace mode).
    pub fn report(self: *TrackingAllocator, writer: anytype) !void {
        const debug_info = try std.debug.getSelfDebugInfo();

        var entries = std.ArrayList(SiteStats).init(self.backing);
        defer entries.deinit();

        var it = self.per_site.valueIterator();
        while (it.next()) |v| try entries.append(v.*);

        std.mem.sort(SiteStats, entries.items, {}, struct {
            fn lt(_: void, a: SiteStats, b: SiteStats) bool {
                return a.peak_live > b.peak_live;
            }
        }.lt);

        const mode = if (self.config.capture_full_trace) "full trace" else "single frame";
        try writer.print(
            "=== Tracking allocator: {d} call sites ({s}) ===\n",
            .{ entries.items.len, mode },
        );

        for (entries.items) |s| {
            if (self.config.capture_full_trace) {
                try writer.print(
                    "\npeak {d} B, live {d} B, alloc'd {d} B, freed {d} B ({d} calls)\n",
                    .{ s.peak_live, s.bytes_live, s.bytes_allocated, s.bytes_freed, s.count },
                );
                const st = std.builtin.StackTrace{
                    .instruction_addresses = @constCast(&s.trace),
                    .index = s.trace_len,
                };
                try std.debug.writeStackTrace(st, writer, debug_info, .no_color);
            } else {
                try writer.print(
                    "peak {d:>10} B  live {d:>10} B  alloc'd {d:>10} B  freed {d:>10} B  ({d:>4} calls)  @ ",
                    .{ s.peak_live, s.bytes_live, s.bytes_allocated, s.bytes_freed, s.count },
                );
                try printAddress(writer, debug_info, s.trace[0]);
                try writer.writeByte('\n');
            }
        }
    }

    fn printAddress(writer: anytype, debug_info: *std.debug.SelfInfo, addr: usize) !void {
        const module = debug_info.getModuleForAddress(addr) catch {
            try writer.print("0x{x}", .{addr});
            return;
        };
        const sym = module.getSymbolAtAddress(debug_info.allocator, addr) catch {
            try writer.print("0x{x}", .{addr});
            return;
        };
        defer if (@hasDecl(@TypeOf(sym), "deinit")) sym.deinit(debug_info.allocator);

        if (sym.source_location) |loc| {
            try writer.print("{s}:{d}:{d} in {s}", .{
                loc.file_name, loc.line, loc.column, sym.name,
            });
        } else {
            try writer.print("{s}", .{sym.name});
        }
    }
};

// --------------- tests -------------------------------------------------------

test "tracks allocations per site (single frame)" {
    const testing = std.testing;

    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();

    var tracker = TrackingAllocator.init(std.heap.page_allocator, gpa.allocator(), .{});
    defer tracker.deinit();

    const a = tracker.allocator();
    const buf1 = try a.alloc(u8, 100);
    const buf2 = try a.alloc(u8, 200);
    defer a.free(buf1);
    defer a.free(buf2);

    try testing.expect(tracker.totalLive() == 300);
}

test "tracks allocations per site (full trace)" {
    const testing = std.testing;

    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();

    var tracker = TrackingAllocator.init(
        std.heap.page_allocator,
        gpa.allocator(),
        .{ .capture_full_trace = true },
    );
    defer tracker.deinit();

    const a = tracker.allocator();
    const buf = try a.alloc(u8, 500);
    try testing.expect(tracker.totalLive() == 500);
    a.free(buf);
    try testing.expect(tracker.totalLive() == 0);
}
