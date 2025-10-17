#pragma once

#include "planar_map.h"

namespace utils {
bool is_terminus(
    const Arrangement_with_history_2 &arr,
    const std::map<Point_2, bool> &point_is_singularity,
    Arrangement_with_history_2::Vertex_const_handle &vertex
);

// assign qi to each segment
std::tuple<bool, std::map<Segment_2, int, Segment2Comparator>, std::vector<Point_2>> validity_check(
  Arrangement_with_history_2 &arr,
  std::map<Point_2, bool> &point_is_singularity,
  // convex/concave labeling per segment. This is based on the "new" segments on the planar map
  std::map<Segment_2, bool, Segment2Comparator> &segment_is_convex,
  // the segment is a cut on the triangulation
  std::map<Segment_2, bool, Segment2Comparator> &segment_is_cut,
  // given an intersection point, which "new" edges are the upper and lower casing edges?
  std::map<Point_2, std::vector<Segment_2>> &upper_casing_edges,
  std::map<Point_2, std::vector<Segment_2>> &lower_casing_edges,
  // given a "new" segment, does segment.source() -> segment.target() matches the canonical orientation?
  std::map<Segment_2, int, Segment2Comparator> &segment_orientation,
  bool is_back_facing
);
} // namespace utils
