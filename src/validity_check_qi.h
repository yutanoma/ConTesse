#pragma once

#include "planar_map.h"

namespace utils {
bool is_terminus(
    const Arrangement_with_history_2 &arr,
    const std::map<Point_2, bool> &point_is_singularity,
    Arrangement_with_history_2::Vertex_const_handle &vertex
);

// assign qi to each segment
std::tuple<bool, std::map<Segment_2, int, Segment2Comparator>, std::vector<Point_2>> validity_check_assign_qi(
    Arrangement_with_history_2 &arr,
    std::map<Point_2, bool> &point_is_singularity,
    std::map<Segment_2, bool, Segment2Comparator> &segment_is_convex,
    // given an intersection point, which edges are the upper and lower casing edges?
    std::map<Point_2, Segment_2> &upper_casing_edges,
    std::map<Point_2, Segment_2> &lower_casing_edges,
    // given a point, which point is the next point in the circular order?
    // *exclude* the self-intersection point in this list
    std::map<Point_2, Point_2> &next
);
} // namespace utils
