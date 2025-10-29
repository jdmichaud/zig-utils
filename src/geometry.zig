const std = @import("std");

/// Reorder vertices so that they form a simple polygon (no intersecting edges)
pub fn make_polygon(vertices: [][2]f32) void {
  const CmpContext = struct { centroid: [2]f32 };

  const Fns = struct {
    // compute in which quadrant the vertex belongs to w.r.t. the centroid
    fn quadrant(vertex: [2]f32, centroid: [2]f32) u8 {
      const dx = vertex[0] - centroid[0];
      const dy = vertex[1] - centroid[1];
      if (dx < 0 and dy < 0) {
        return 0;
      } else if (dx < 0 and dy >= 0) {
        return 1;
      } else if (dx >= 0 and dy >= 0) {
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
      const dx = vertex[0] - centroid[0];
      const dy = vertex[1] - centroid[1];
      return dx * dx + dy * dy;
    }

    // order by quadrant, then by cross product and then by distance from centroid
    fn cmpfn(context: CmpContext, a: [2]f32, b: [2]f32) bool {
      const qa = quadrant(a, context.centroid);
      const qb = quadrant(b, context.centroid);
      if (qa != qb) {
        return qa < qb;
      }
      // in same quadrant: sort by cross product
      const c = cross(a, b, context.centroid);
      if (c != 0) {
        return c < 0;
      }
      // Collinear points: sort by distance from centroid
      const d1 = dist_sq(a, context.centroid);
      const d2 = dist_sq(b, context.centroid);
      return d1 > d2;
    }
  };

  var average_x: f32 = 0;
  var average_y: f32 = 0;
  for (vertices) |vertex| {
    average_x += vertex[0];
    average_y += vertex[1];
  }
  average_x /= @floatFromInt(vertices.len) ;
  average_y /= @floatFromInt(vertices.len);
  const centroid: [2]f32 = .{ average_x, average_y };

  std.mem.sort([2]f32, vertices, CmpContext{ .centroid = centroid }, Fns.cmpfn);
}

test "make_polygon" {
  std.testing.log_level = .debug;

  {
    var vertices = [_][2]f32{ [2]f32{ 0.0, 0.0 }, [2]f32{ 1.0, 1.0 }, [2]f32{ 1.0, 0.0 }, [2]f32{ 0.0, 1.0 } };
    make_polygon(&vertices);
    try std.testing.expectEqual(.{ 0.0, 0.0 }, vertices[0]);
    try std.testing.expectEqual(.{ 0.0, 1.0 }, vertices[1]);
    try std.testing.expectEqual(.{ 1.0, 1.0 }, vertices[2]);
    try std.testing.expectEqual(.{ 1.0, 0.0 }, vertices[3]);
  }
  {
    var vertices = [_][2]f32{ [2]f32{ 0.0, 0.0 }, [2]f32{ 1.0, 0.0 }, [2]f32{ 0.0, 1.0 }, [2]f32{ 5.0, 5.0 } };
    make_polygon(&vertices);
    try std.testing.expectEqual(.{ 1.0, 0.0 }, vertices[0]);
    try std.testing.expectEqual(.{ 0.0, 0.0 }, vertices[1]);
    try std.testing.expectEqual(.{ 0.0, 1.0 }, vertices[2]);
    try std.testing.expectEqual(.{ 5.0, 5.0 }, vertices[3]);
  }
}

