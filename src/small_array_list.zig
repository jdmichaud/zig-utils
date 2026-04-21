//! A contiguous, growable list of items that uses a caller-provided static
//! buffer for small-buffer optimization. The list writes into the static
//! buffer until it is full, then falls back to allocating from the given
//! `Allocator`. Once on the heap, `shrinkAndFree` may move items back into
//! the static buffer if they fit.
//!
//! The public API mirrors `std.ArrayList`. `init` takes a single allocator
//! (just like `std.ArrayList`); use `initFromBuffer` to provide a static
//! buffer for small-buffer optimization.
//!
//! Written with Claude Opus 4.7.

const std = @import("std");
const assert = std.debug.assert;
const testing = std.testing;
const mem = std.mem;
const math = std.math;
const Allocator = mem.Allocator;

pub fn SmallArrayList(comptime T: type) type {
  return struct {
    const Self = @This();

    /// Contents of the list. May be accessed directly.
    items: Slice,
    /// How many T values this list can hold without allocating more memory.
    capacity: usize,
    allocator: Allocator,
    /// Externally-provided static buffer. When `items.ptr == static_buffer.ptr`
    /// the list is using the static buffer and must not be freed.
    static_buffer: Slice,

    pub const Slice = []T;

    pub fn SentinelSlice(comptime s: T) type {
      return [:s]T;
    }

    /// Initialize with a caller-provided static buffer. The list will use
    /// this buffer until it grows beyond `buffer.len`, at which point it
    /// allocates heap memory. Passing an empty slice disables SBO.
    pub fn initFromBuffer(allocator: Allocator, buffer: Slice) Self {
      return .{
        .items = buffer[0..0],
        .capacity = buffer.len,
        .allocator = allocator,
        .static_buffer = buffer,
      };
    }

    /// Initialize without a static buffer (behaves like a regular ArrayList).
    pub fn init(allocator: Allocator) Self {
      return .{
        .items = &[_]T{},
        .capacity = 0,
        .allocator = allocator,
        .static_buffer = &[_]T{},
      };
    }

    /// Initialize with capacity to hold `num` elements. No static buffer.
    pub fn initCapacity(allocator: Allocator, num: usize) Allocator.Error!Self {
      var self = Self.init(allocator);
      try self.ensureTotalCapacityPrecise(num);
      return self;
    }

    /// True if the current items live inside the static buffer.
    pub fn isStatic(self: Self) bool {
      if (self.static_buffer.len == 0) return false;
      return @intFromPtr(self.items.ptr) == @intFromPtr(self.static_buffer.ptr);
    }

    /// Release all allocated memory (no-op if using the static buffer).
    pub fn deinit(self: Self) void {
      if (@sizeOf(T) > 0 and !self.isStatic() and self.capacity > 0) {
        self.allocator.free(self.allocatedSlice());
      }
    }

    /// Takes ownership of the passed in slice. The slice must have been
    /// allocated with `allocator`. No static buffer.
    pub fn fromOwnedSlice(allocator: Allocator, slice: Slice) Self {
      return .{
        .items = slice,
        .capacity = slice.len,
        .allocator = allocator,
        .static_buffer = &[_]T{},
      };
    }

    pub fn fromOwnedSliceSentinel(
      allocator: Allocator,
      comptime sentinel: T,
      slice: [:sentinel]T,
    ) Self {
      return .{
        .items = slice[0 .. slice.len + 1],
        .capacity = slice.len + 1,
        .allocator = allocator,
        .static_buffer = &[_]T{},
      };
    }

    /// The caller owns the returned memory. Empties this list.
    /// If the items currently live in the static buffer, a new heap
    /// allocation is made and the contents are copied.
    pub fn toOwnedSlice(self: *Self) Allocator.Error!Slice {
      const allocator = self.allocator;
      const static_buf = self.static_buffer;

      if (self.isStatic()) {
        const new_memory = try allocator.alloc(T, self.items.len);
        @memcpy(new_memory, self.items);
        self.* = Self.initFromBuffer(allocator, static_buf);
        return new_memory;
      }

      const old_memory = self.allocatedSlice();
      if (allocator.remap(old_memory, self.items.len)) |new_items| {
        self.* = Self.initFromBuffer(allocator, static_buf);
        return new_items;
      }

      const new_memory = try allocator.alloc(T, self.items.len);
      @memcpy(new_memory, self.items);
      self.clearAndFree();
      return new_memory;
    }

    pub fn toOwnedSliceSentinel(
      self: *Self,
      comptime sentinel: T,
    ) Allocator.Error!SentinelSlice(sentinel) {
      try self.ensureTotalCapacityPrecise(self.items.len + 1);
      self.appendAssumeCapacity(sentinel);
      const result = try self.toOwnedSlice();
      return result[0 .. result.len - 1 :sentinel];
    }

    /// Creates a copy of this list. The clone is always heap-backed
    /// (never reuses the original's static buffer).
    pub fn clone(self: Self) Allocator.Error!Self {
      var cloned = Self.init(self.allocator);
      errdefer cloned.deinit();
      try cloned.ensureTotalCapacityPrecise(self.capacity);
      cloned.appendSliceAssumeCapacity(self.items);
      return cloned;
    }

    pub fn insert(self: *Self, i: usize, item: T) Allocator.Error!void {
      const dst = try self.addManyAt(i, 1);
      dst[0] = item;
    }

    pub fn insertAssumeCapacity(self: *Self, i: usize, item: T) void {
      assert(self.items.len < self.capacity);
      self.items.len += 1;
      mem.copyBackwards(T, self.items[i + 1 .. self.items.len], self.items[i .. self.items.len - 1]);
      self.items[i] = item;
    }

    pub fn addManyAt(self: *Self, index: usize, count: usize) Allocator.Error![]T {
      const new_len = try addOrOom(self.items.len, count);
      if (self.capacity >= new_len)
        return addManyAtAssumeCapacity(self, index, count);

      const new_capacity = growCapacity(self.capacity, new_len);

      // Try to grow in place (only possible if we're on the heap).
      if (!self.isStatic() and self.capacity > 0) {
        const old_memory = self.allocatedSlice();
        if (self.allocator.remap(old_memory, new_capacity)) |new_memory| {
          self.items.ptr = new_memory.ptr;
          self.capacity = new_memory.len;
          return addManyAtAssumeCapacity(self, index, count);
        }
      }

      const new_memory = try self.allocator.alloc(T, new_capacity);
      const to_move = self.items[index..];
      @memcpy(new_memory[0..index], self.items[0..index]);
      @memcpy(new_memory[index + count ..][0..to_move.len], to_move);

      // Free previous heap allocation (if any). Never free the static buffer.
      if (!self.isStatic() and self.capacity > 0) {
        self.allocator.free(self.allocatedSlice());
      }

      self.items = new_memory[0..new_len];
      self.capacity = new_memory.len;
      return new_memory[index..][0..count];
    }

    pub fn addManyAtAssumeCapacity(self: *Self, index: usize, count: usize) []T {
      const new_len = self.items.len + count;
      assert(self.capacity >= new_len);
      const to_move = self.items[index..];
      self.items.len = new_len;
      mem.copyBackwards(T, self.items[index + count ..], to_move);
      const result = self.items[index..][0..count];
      @memset(result, undefined);
      return result;
    }

    pub fn insertSlice(self: *Self, index: usize, items: []const T) Allocator.Error!void {
      const dst = try self.addManyAt(index, items.len);
      @memcpy(dst, items);
    }

    pub fn replaceRange(self: *Self, start: usize, len: usize, new_items: []const T) Allocator.Error!void {
      const after_range = start + len;
      const range = self.items[start..after_range];
      if (range.len < new_items.len) {
        const first = new_items[0..range.len];
        const rest = new_items[range.len..];
        @memcpy(range[0..first.len], first);
        try self.insertSlice(after_range, rest);
      } else {
        self.replaceRangeAssumeCapacity(start, len, new_items);
      }
    }

    pub fn replaceRangeAssumeCapacity(self: *Self, start: usize, len: usize, new_items: []const T) void {
      const after_range = start + len;
      const range = self.items[start..after_range];

      if (range.len == new_items.len) {
        @memcpy(range[0..new_items.len], new_items);
      } else if (range.len < new_items.len) {
        const first = new_items[0..range.len];
        const rest = new_items[range.len..];
        @memcpy(range[0..first.len], first);
        const dst = self.addManyAtAssumeCapacity(after_range, rest.len);
        @memcpy(dst, rest);
      } else {
        const extra = range.len - new_items.len;
        @memcpy(range[0..new_items.len], new_items);
        std.mem.copyForwards(
          T,
          self.items[after_range - extra ..],
          self.items[after_range..],
        );
        @memset(self.items[self.items.len - extra ..], undefined);
        self.items.len -= extra;
      }
    }

    pub fn append(self: *Self, item: T) Allocator.Error!void {
      const new_item_ptr = try self.addOne();
      new_item_ptr.* = item;
    }

    pub fn appendAssumeCapacity(self: *Self, item: T) void {
      self.addOneAssumeCapacity().* = item;
    }

    pub fn orderedRemove(self: *Self, i: usize) T {
      const old_item = self.items[i];
      self.replaceRangeAssumeCapacity(i, 1, &.{});
      return old_item;
    }

    pub fn swapRemove(self: *Self, i: usize) T {
      if (self.items.len - 1 == i) return self.pop().?;
      const old_item = self.items[i];
      self.items[i] = self.pop().?;
      return old_item;
    }

    pub fn appendSlice(self: *Self, items_to_add: []const T) Allocator.Error!void {
      try self.ensureUnusedCapacity(items_to_add.len);
      self.appendSliceAssumeCapacity(items_to_add);
    }

    pub fn appendSliceAssumeCapacity(self: *Self, items_to_add: []const T) void {
      const old_len = self.items.len;
      const new_len = old_len + items_to_add.len;
      assert(new_len <= self.capacity);
      self.items.len = new_len;
      @memcpy(self.items[old_len..][0..items_to_add.len], items_to_add);
    }

    pub fn appendUnalignedSlice(self: *Self, items_to_add: []align(1) const T) Allocator.Error!void {
      try self.ensureUnusedCapacity(items_to_add.len);
      self.appendUnalignedSliceAssumeCapacity(items_to_add);
    }

    pub fn appendUnalignedSliceAssumeCapacity(self: *Self, items_to_add: []align(1) const T) void {
      const old_len = self.items.len;
      const new_len = old_len + items_to_add.len;
      assert(new_len <= self.capacity);
      self.items.len = new_len;
      @memcpy(self.items[old_len..][0..items_to_add.len], items_to_add);
    }

    pub inline fn appendNTimes(self: *Self, value: T, n: usize) Allocator.Error!void {
      const old_len = self.items.len;
      try self.resize(try addOrOom(old_len, n));
      @memset(self.items[old_len..self.items.len], value);
    }

    pub inline fn appendNTimesAssumeCapacity(self: *Self, value: T, n: usize) void {
      const new_len = self.items.len + n;
      assert(new_len <= self.capacity);
      @memset(self.items.ptr[self.items.len..new_len], value);
      self.items.len = new_len;
    }

    pub fn resize(self: *Self, new_len: usize) Allocator.Error!void {
      try self.ensureTotalCapacity(new_len);
      self.items.len = new_len;
    }

    pub fn shrinkAndFree(self: *Self, new_len: usize) void {
      assert(new_len <= self.items.len);

      if (@sizeOf(T) == 0) {
        self.items.len = new_len;
        return;
      }

      if (self.isStatic()) {
        // Already in the static buffer; just adjust length.
        self.items.len = new_len;
        return;
      }

      if (self.capacity == 0) {
        self.items.len = new_len;
        return;
      }

      // If new_len fits in the static buffer, migrate back and free heap.
      if (new_len <= self.static_buffer.len) {
        @memcpy(self.static_buffer[0..new_len], self.items[0..new_len]);
        const old_memory = self.allocatedSlice();
        self.items = self.static_buffer[0..new_len];
        self.capacity = self.static_buffer.len;
        self.allocator.free(old_memory);
        return;
      }

      const old_memory = self.allocatedSlice();
      if (self.allocator.remap(old_memory, new_len)) |new_items| {
        self.capacity = new_items.len;
        self.items.ptr = new_items.ptr;
        self.items.len = new_len;
        return;
      }

      const new_memory = self.allocator.alloc(T, new_len) catch |e| switch (e) {
        error.OutOfMemory => {
          self.items.len = new_len;
          return;
        },
      };
      @memcpy(new_memory, self.items[0..new_len]);
      self.allocator.free(old_memory);
      self.items = new_memory;
      self.capacity = new_memory.len;
    }

    pub fn shrinkRetainingCapacity(self: *Self, new_len: usize) void {
      assert(new_len <= self.items.len);
      self.items.len = new_len;
    }

    pub fn clearRetainingCapacity(self: *Self) void {
      self.items.len = 0;
    }

    pub fn clearAndFree(self: *Self) void {
      if (!self.isStatic() and self.capacity > 0) {
        self.allocator.free(self.allocatedSlice());
      }
      if (self.static_buffer.len > 0) {
        self.items = self.static_buffer[0..0];
        self.capacity = self.static_buffer.len;
      } else {
        self.items = &[_]T{};
        self.capacity = 0;
      }
    }

    pub fn ensureTotalCapacity(self: *Self, new_capacity: usize) Allocator.Error!void {
      if (@sizeOf(T) == 0) {
        self.capacity = math.maxInt(usize);
        return;
      }
      if (self.capacity >= new_capacity) return;
      const better_capacity = growCapacity(self.capacity, new_capacity);
      return self.ensureTotalCapacityPrecise(better_capacity);
    }

    pub fn ensureTotalCapacityPrecise(self: *Self, new_capacity: usize) Allocator.Error!void {
      if (@sizeOf(T) == 0) {
        self.capacity = math.maxInt(usize);
        return;
      }
      if (self.capacity >= new_capacity) return;

      if (self.isStatic()) {
        // Migrate from static buffer to heap.
        const new_memory = try self.allocator.alloc(T, new_capacity);
        @memcpy(new_memory[0..self.items.len], self.items);
        self.items.ptr = new_memory.ptr;
        self.capacity = new_memory.len;
        return;
      }

      if (self.capacity == 0) {
        const new_memory = try self.allocator.alloc(T, new_capacity);
        self.items.ptr = new_memory.ptr;
        self.capacity = new_memory.len;
        return;
      }

      const old_memory = self.allocatedSlice();
      if (self.allocator.remap(old_memory, new_capacity)) |new_memory| {
        self.items.ptr = new_memory.ptr;
        self.capacity = new_memory.len;
      } else {
        const new_memory = try self.allocator.alloc(T, new_capacity);
        @memcpy(new_memory[0..self.items.len], self.items);
        self.allocator.free(old_memory);
        self.items.ptr = new_memory.ptr;
        self.capacity = new_memory.len;
      }
    }

    pub fn ensureUnusedCapacity(self: *Self, additional_count: usize) Allocator.Error!void {
      return self.ensureTotalCapacity(try addOrOom(self.items.len, additional_count));
    }

    pub fn expandToCapacity(self: *Self) void {
      self.items.len = self.capacity;
    }

    pub fn addOne(self: *Self) Allocator.Error!*T {
      const newlen = self.items.len + 1;
      try self.ensureTotalCapacity(newlen);
      return self.addOneAssumeCapacity();
    }

    pub fn addOneAssumeCapacity(self: *Self) *T {
      assert(self.items.len < self.capacity);
      self.items.len += 1;
      return &self.items[self.items.len - 1];
    }

    pub fn addManyAsArray(self: *Self, comptime n: usize) Allocator.Error!*[n]T {
      const prev_len = self.items.len;
      try self.resize(try addOrOom(self.items.len, n));
      return self.items[prev_len..][0..n];
    }

    pub fn addManyAsArrayAssumeCapacity(self: *Self, comptime n: usize) *[n]T {
      assert(self.items.len + n <= self.capacity);
      const prev_len = self.items.len;
      self.items.len += n;
      return self.items[prev_len..][0..n];
    }

    pub fn addManyAsSlice(self: *Self, n: usize) Allocator.Error![]T {
      const prev_len = self.items.len;
      try self.resize(try addOrOom(self.items.len, n));
      return self.items[prev_len..][0..n];
    }

    pub fn addManyAsSliceAssumeCapacity(self: *Self, n: usize) []T {
      assert(self.items.len + n <= self.capacity);
      const prev_len = self.items.len;
      self.items.len += n;
      return self.items[prev_len..][0..n];
    }

    pub fn pop(self: *Self) ?T {
      if (self.items.len == 0) return null;
      const val = self.items[self.items.len - 1];
      self.items.len -= 1;
      return val;
    }

    pub fn allocatedSlice(self: Self) Slice {
      return self.items.ptr[0..self.capacity];
    }

    pub fn unusedCapacitySlice(self: Self) []T {
      return self.allocatedSlice()[self.items.len..];
    }

    pub fn getLast(self: Self) T {
      return self.items[self.items.len - 1];
    }

    pub fn getLastOrNull(self: Self) ?T {
      if (self.items.len == 0) return null;
      return self.getLast();
    }

    const init_capacity: usize = @max(1, std.atomic.cache_line / @sizeOf(T));

    fn growCapacity(current: usize, minimum: usize) usize {
      var new = current;
      while (true) {
        new +|= new / 2 + init_capacity;
        if (new >= minimum) return new;
      }
    }
  };
}

fn addOrOom(a: usize, b: usize) error{OutOfMemory}!usize {
  const result, const overflow = @addWithOverflow(a, b);
  if (overflow != 0) return error.OutOfMemory;
  return result;
}

// -----------------------------------------------------------------------------
// Tests: these mirror the tests in array_list.zig, adapted for SmallArrayList.
// Each test runs twice (where useful): once with a static buffer large enough
// to hold everything (exercising the static path) and once with a small or
// empty buffer (exercising the heap fallback and the migration between them).
// -----------------------------------------------------------------------------

test "init" {
  var buf: [4]i32 = undefined;
  var list = SmallArrayList(i32).initFromBuffer(testing.allocator, &buf);
  defer list.deinit();

  try testing.expect(list.items.len == 0);
  try testing.expect(list.capacity == 4);
  try testing.expect(list.isStatic());
}

test "init without static buffer" {
  var list = SmallArrayList(i32).init(testing.allocator);
  defer list.deinit();

  try testing.expect(list.items.len == 0);
  try testing.expect(list.capacity == 0);
  try testing.expect(!list.isStatic());
}

test "initCapacity" {
  const a = testing.allocator;
  var list = try SmallArrayList(i8).initCapacity(a, 200);
  defer list.deinit();
  try testing.expect(list.items.len == 0);
  try testing.expect(list.capacity >= 200);
  try testing.expect(!list.isStatic());
}

test "clone" {
  const a = testing.allocator;
  var buf: [8]i32 = undefined;
  var array = SmallArrayList(i32).initFromBuffer(a, &buf);
  try array.append(-1);
  try array.append(3);
  try array.append(5);

  var cloned = try array.clone();
  defer cloned.deinit();

  try testing.expectEqualSlices(i32, array.items, cloned.items);
  try testing.expectEqual(array.allocator.ptr, cloned.allocator.ptr);
  try testing.expect(cloned.capacity >= array.capacity);

  array.deinit();

  try testing.expectEqual(@as(i32, -1), cloned.items[0]);
  try testing.expectEqual(@as(i32, 3), cloned.items[1]);
  try testing.expectEqual(@as(i32, 5), cloned.items[2]);
}

fn runBasic(buffer_size: comptime_int) !void {
  const a = testing.allocator;
  var buf: [buffer_size]i32 = undefined;
  var list = SmallArrayList(i32).initFromBuffer(a, &buf);
  defer list.deinit();

  {
    var i: usize = 0;
    while (i < 10) : (i += 1) {
      list.append(@as(i32, @intCast(i + 1))) catch unreachable;
    }
  }

  {
    var i: usize = 0;
    while (i < 10) : (i += 1) {
      try testing.expect(list.items[i] == @as(i32, @intCast(i + 1)));
    }
  }

  for (list.items, 0..) |v, i| {
    try testing.expect(v == @as(i32, @intCast(i + 1)));
  }

  try testing.expect(list.pop() == 10);
  try testing.expect(list.items.len == 9);

  list.appendSlice(&[_]i32{ 1, 2, 3 }) catch unreachable;
  try testing.expect(list.items.len == 12);
  try testing.expect(list.pop() == 3);
  try testing.expect(list.pop() == 2);
  try testing.expect(list.pop() == 1);
  try testing.expect(list.items.len == 9);

  var unaligned: [3]i32 align(1) = [_]i32{ 4, 5, 6 };
  list.appendUnalignedSlice(&unaligned) catch unreachable;
  try testing.expect(list.items.len == 12);
  try testing.expect(list.pop() == 6);
  try testing.expect(list.pop() == 5);
  try testing.expect(list.pop() == 4);
  try testing.expect(list.items.len == 9);

  list.appendSlice(&[_]i32{}) catch unreachable;
  try testing.expect(list.items.len == 9);

  list.items[7] = 33;
  list.items[8] = 42;

  try testing.expect(list.pop() == 42);
  try testing.expect(list.pop() == 33);
}

test "basic (stays in static buffer)" {
  try runBasic(32);
}

test "basic (overflows to heap)" {
  try runBasic(4);
}

test "basic (no static buffer)" {
  try runBasic(0);
}

test "appendNTimes" {
  const a = testing.allocator;
  {
    var buf: [16]i32 = undefined;
    var list = SmallArrayList(i32).initFromBuffer(a, &buf);
    defer list.deinit();

    try list.appendNTimes(2, 10);
    try testing.expectEqual(@as(usize, 10), list.items.len);
    try testing.expect(list.isStatic());
    for (list.items) |element| {
      try testing.expectEqual(@as(i32, 2), element);
    }
  }
  {
    var buf: [4]i32 = undefined;
    var list = SmallArrayList(i32).initFromBuffer(a, &buf);
    defer list.deinit();

    try list.appendNTimes(2, 10);
    try testing.expectEqual(@as(usize, 10), list.items.len);
    try testing.expect(!list.isStatic());
    for (list.items) |element| {
      try testing.expectEqual(@as(i32, 2), element);
    }
  }
}

test "appendNTimes with failing allocator" {
  const a = testing.failing_allocator;
  var list = SmallArrayList(i32).init(a);
  defer list.deinit();
  try testing.expectError(error.OutOfMemory, list.appendNTimes(2, 10));
}

test "orderedRemove" {
  const a = testing.allocator;
  inline for (.{ 32, 4 }) |N| {
    var buf: [N]i32 = undefined;
    var list = SmallArrayList(i32).initFromBuffer(a, &buf);
    defer list.deinit();

    try list.append(1);
    try list.append(2);
    try list.append(3);
    try list.append(4);
    try list.append(5);
    try list.append(6);
    try list.append(7);

    try testing.expectEqual(@as(i32, 4), list.orderedRemove(3));
    try testing.expectEqual(@as(i32, 5), list.items[3]);
    try testing.expectEqual(@as(usize, 6), list.items.len);

    try testing.expectEqual(@as(i32, 7), list.orderedRemove(5));
    try testing.expectEqual(@as(usize, 5), list.items.len);

    try testing.expectEqual(@as(i32, 1), list.orderedRemove(0));
    try testing.expectEqual(@as(i32, 2), list.items[0]);
    try testing.expectEqual(@as(usize, 4), list.items.len);
  }

  // remove last item
  {
    var buf: [4]i32 = undefined;
    var list = SmallArrayList(i32).initFromBuffer(a, &buf);
    defer list.deinit();
    try list.append(1);
    try testing.expectEqual(@as(i32, 1), list.orderedRemove(0));
    try testing.expectEqual(@as(usize, 0), list.items.len);
  }
}

test "swapRemove" {
  const a = testing.allocator;
  inline for (.{ 32, 4 }) |N| {
    var buf: [N]i32 = undefined;
    var list = SmallArrayList(i32).initFromBuffer(a, &buf);
    defer list.deinit();

    try list.append(1);
    try list.append(2);
    try list.append(3);
    try list.append(4);
    try list.append(5);
    try list.append(6);
    try list.append(7);

    try testing.expect(list.swapRemove(3) == 4);
    try testing.expect(list.items[3] == 7);
    try testing.expect(list.items.len == 6);

    try testing.expect(list.swapRemove(5) == 6);
    try testing.expect(list.items.len == 5);

    try testing.expect(list.swapRemove(0) == 1);
    try testing.expect(list.items[0] == 5);
    try testing.expect(list.items.len == 4);
  }
}

test "insert" {
  const a = testing.allocator;
  inline for (.{ 16, 2 }) |N| {
    var buf: [N]i32 = undefined;
    var list = SmallArrayList(i32).initFromBuffer(a, &buf);
    defer list.deinit();

    try list.insert(0, 1);
    try list.append(2);
    try list.insert(2, 3);
    try list.insert(0, 5);
    try testing.expect(list.items[0] == 5);
    try testing.expect(list.items[1] == 1);
    try testing.expect(list.items[2] == 2);
    try testing.expect(list.items[3] == 3);
  }
}

test "insertSlice" {
  const a = testing.allocator;
  inline for (.{ 16, 4 }) |N| {
    var buf: [N]i32 = undefined;
    var list = SmallArrayList(i32).initFromBuffer(a, &buf);
    defer list.deinit();

    try list.append(1);
    try list.append(2);
    try list.append(3);
    try list.append(4);
    try list.insertSlice(1, &[_]i32{ 9, 8 });
    try testing.expect(list.items[0] == 1);
    try testing.expect(list.items[1] == 9);
    try testing.expect(list.items[2] == 8);
    try testing.expect(list.items[3] == 2);
    try testing.expect(list.items[4] == 3);
    try testing.expect(list.items[5] == 4);

    const items = [_]i32{1};
    try list.insertSlice(0, items[0..0]);
    try testing.expect(list.items.len == 6);
    try testing.expect(list.items[0] == 1);
  }
}

test "replaceRange" {
  const a = testing.allocator;
  const Case = struct {
    start: usize,
    len: usize,
    new_items: []const i32,
    expected: []const i32,
  };
  const cases = [_]Case{
    .{ .start = 1, .len = 0, .new_items = &.{ 0, 0, 0 }, .expected = &.{ 1, 0, 0, 0, 2, 3, 4, 5 } },
    .{ .start = 1, .len = 1, .new_items = &.{ 0, 0, 0 }, .expected = &.{ 1, 0, 0, 0, 3, 4, 5 } },
    .{ .start = 1, .len = 2, .new_items = &.{ 0, 0, 0 }, .expected = &.{ 1, 0, 0, 0, 4, 5 } },
    .{ .start = 1, .len = 3, .new_items = &.{ 0, 0, 0 }, .expected = &.{ 1, 0, 0, 0, 5 } },
    .{ .start = 1, .len = 4, .new_items = &.{ 0, 0, 0 }, .expected = &.{ 1, 0, 0, 0 } },
  };

  inline for (.{ 16, 4 }) |N| {
    for (cases) |c| {
      var buf: [N]i32 = undefined;
      var list = SmallArrayList(i32).initFromBuffer(a, &buf);
      defer list.deinit();
      try list.appendSlice(&[_]i32{ 1, 2, 3, 4, 5 });
      try list.replaceRange(c.start, c.len, c.new_items);
      try testing.expectEqualSlices(i32, c.expected, list.items);
    }
  }
}

test "replaceRangeAssumeCapacity" {
  const a = testing.allocator;
  const cases = .{
    .{ .start = 1, .len = 0, .new_items = &[_]i32{ 0, 0, 0 }, .expected = &[_]i32{ 1, 0, 0, 0, 2, 3, 4, 5 } },
    .{ .start = 1, .len = 1, .new_items = &[_]i32{ 0, 0, 0 }, .expected = &[_]i32{ 1, 0, 0, 0, 3, 4, 5 } },
    .{ .start = 1, .len = 2, .new_items = &[_]i32{ 0, 0, 0 }, .expected = &[_]i32{ 1, 0, 0, 0, 4, 5 } },
    .{ .start = 1, .len = 3, .new_items = &[_]i32{ 0, 0, 0 }, .expected = &[_]i32{ 1, 0, 0, 0, 5 } },
    .{ .start = 1, .len = 4, .new_items = &[_]i32{ 0, 0, 0 }, .expected = &[_]i32{ 1, 0, 0, 0 } },
  };

  inline for (cases) |c| {
    var buf: [16]i32 = undefined;
    var list = SmallArrayList(i32).initFromBuffer(a, &buf);
    defer list.deinit();
    try list.appendSlice(&[_]i32{ 1, 2, 3, 4, 5 });
    list.replaceRangeAssumeCapacity(c.start, c.len, c.new_items);
    try testing.expectEqualSlices(i32, c.expected, list.items);
  }
}

test "SmallArrayList(T) of struct T" {
  const a = std.testing.allocator;
  const Item = struct {
    integer: i32,
    sub_items: SmallArrayList(@This()),
  };
  var root = Item{
    .integer = 1,
    .sub_items = SmallArrayList(Item).init(a),
  };
  defer root.sub_items.deinit();

  try root.sub_items.append(Item{
    .integer = 42,
    .sub_items = SmallArrayList(Item).init(a),
  });
  try testing.expect(root.sub_items.items[0].integer == 42);
}

test "shrink still sets length when resizing is disabled" {
  var failing_allocator = testing.FailingAllocator.init(testing.allocator, .{ .resize_fail_index = 0 });
  const a = failing_allocator.allocator();

  // On heap — resize may fail but shrink must still adjust length.
  var list = SmallArrayList(i32).init(a);
  defer list.deinit();

  try list.append(1);
  try list.append(2);
  try list.append(3);

  list.shrinkAndFree(1);
  try testing.expect(list.items.len == 1);
}

test "shrinkAndFree with a copy" {
  var failing_allocator = testing.FailingAllocator.init(testing.allocator, .{ .resize_fail_index = 0 });
  const a = failing_allocator.allocator();

  var list = SmallArrayList(i32).init(a);
  defer list.deinit();

  try list.appendNTimes(3, 16);
  list.shrinkAndFree(4);
  try testing.expect(mem.eql(i32, list.items, &.{ 3, 3, 3, 3 }));
}

test "shrinkAndFree migrates back to static buffer" {
  const a = testing.allocator;
  var buf: [4]i32 = undefined;
  var list = SmallArrayList(i32).initFromBuffer(a, &buf);
  defer list.deinit();

  try list.appendNTimes(7, 16);
  try testing.expect(!list.isStatic());

  list.shrinkAndFree(3);
  try testing.expect(list.isStatic());
  try testing.expectEqualSlices(i32, &.{ 7, 7, 7 }, list.items);
}

test "addManyAsArray" {
  const a = std.testing.allocator;
  inline for (.{ 16, 2 }) |N| {
    var buf: [N]u8 = undefined;
    var list = SmallArrayList(u8).initFromBuffer(a, &buf);
    defer list.deinit();

    (try list.addManyAsArray(4)).* = "aoeu".*;
    try list.ensureTotalCapacity(8);
    list.addManyAsArrayAssumeCapacity(4).* = "asdf".*;

    try testing.expectEqualSlices(u8, list.items, "aoeuasdf");
  }
}

test "growing memory preserves contents" {
  const a = std.testing.allocator;
  inline for (.{ 16, 4, 2 }) |N| {
    var buf: [N]u8 = undefined;
    var list = SmallArrayList(u8).initFromBuffer(a, &buf);
    defer list.deinit();

    (try list.addManyAsArray(4)).* = "abcd".*;
    list.shrinkAndFree(4);

    try list.appendSlice("efgh");
    try testing.expectEqualSlices(u8, list.items, "abcdefgh");
    list.shrinkAndFree(8);

    try list.insertSlice(4, "ijkl");
    try testing.expectEqualSlices(u8, list.items, "abcdijklefgh");
  }
}

test "fromOwnedSlice" {
  const a = testing.allocator;
  var buf: [4]u8 = undefined;
  var orig_list = SmallArrayList(u8).initFromBuffer(a, &buf);
  defer orig_list.deinit();
  try orig_list.appendSlice("foobar");

  const slice = try orig_list.toOwnedSlice();
  var list = SmallArrayList(u8).fromOwnedSlice(a, slice);
  defer list.deinit();
  try testing.expectEqualStrings(list.items, "foobar");
}

test "toOwnedSliceSentinel" {
  const a = testing.allocator;
  inline for (.{ 16, 4 }) |N| {
    var buf: [N]u8 = undefined;
    var list = SmallArrayList(u8).initFromBuffer(a, &buf);
    defer list.deinit();

    try list.appendSlice("foobar");

    const result = try list.toOwnedSliceSentinel(0);
    defer a.free(result);
    try testing.expectEqualStrings(result, mem.sliceTo(result.ptr, 0));
  }
}

test "SmallArrayList(u0)" {
  const a = testing.failing_allocator;

  var list = SmallArrayList(u0).init(a);
  defer list.deinit();

  try list.append(0);
  try list.append(0);
  try list.append(0);
  try testing.expectEqual(list.items.len, 3);

  var count: usize = 0;
  for (list.items) |x| {
    try testing.expectEqual(x, 0);
    count += 1;
  }
  try testing.expectEqual(count, 3);
}

test "SmallArrayList(?u32).pop()" {
  const a = testing.allocator;
  var buf: [4]?u32 = undefined;
  var list = SmallArrayList(?u32).initFromBuffer(a, &buf);
  defer list.deinit();

  try list.append(null);
  try list.append(1);
  try list.append(2);
  try testing.expectEqual(list.items.len, 3);

  try testing.expect(list.pop().? == @as(u32, 2));
  try testing.expect(list.pop().? == @as(u32, 1));
  try testing.expect(list.pop().? == null);
  try testing.expect(list.pop() == null);
}

test "SmallArrayList(u32).getLast()" {
  const a = testing.allocator;
  var buf: [4]u32 = undefined;
  var list = SmallArrayList(u32).initFromBuffer(a, &buf);
  defer list.deinit();

  try list.append(2);
  const const_list = list;
  try testing.expectEqual(const_list.getLast(), 2);
}

test "SmallArrayList(u32).getLastOrNull()" {
  const a = testing.allocator;
  var buf: [4]u32 = undefined;
  var list = SmallArrayList(u32).initFromBuffer(a, &buf);
  defer list.deinit();

  try testing.expectEqual(list.getLastOrNull(), null);

  try list.append(2);
  const const_list = list;
  try testing.expectEqual(const_list.getLastOrNull().?, 2);
}

test "return OutOfMemory when capacity would exceed maximum usize integer value" {
  const a = testing.allocator;
  const new_item: u32 = 42;
  const items = &.{ 42, 43 };

  var list: SmallArrayList(u32) = .{
    .items = undefined,
    .capacity = math.maxInt(usize) - 1,
    .allocator = a,
    .static_buffer = &[_]u32{},
  };
  list.items.len = math.maxInt(usize) - 1;

  try testing.expectError(error.OutOfMemory, list.appendSlice(items));
  try testing.expectError(error.OutOfMemory, list.appendNTimes(new_item, 2));
  try testing.expectError(error.OutOfMemory, list.appendUnalignedSlice(&.{ new_item, new_item }));
  try testing.expectError(error.OutOfMemory, list.addManyAt(0, 2));
  try testing.expectError(error.OutOfMemory, list.addManyAsArray(2));
  try testing.expectError(error.OutOfMemory, list.addManyAsSlice(2));
  try testing.expectError(error.OutOfMemory, list.insertSlice(0, items));
  try testing.expectError(error.OutOfMemory, list.ensureUnusedCapacity(2));
}
