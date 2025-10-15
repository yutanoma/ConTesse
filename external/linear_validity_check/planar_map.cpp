#include "planar_map.h"
#include <map>
#include <iostream>
#include <cmath>
#include <algorithm>

#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_with_history_2.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace utils {
bool is_identical_segment(const Segment_2 &s1, const Segment_2 &s2) {
  return (s1.source() == s2.source() && s1.target() == s2.target()) ||
          (s1.source() == s2.target() && s1.target() == s2.source());
}
  
std::tuple<
    Arrangement_with_history_2,
    std::map<Point_2, int>,
    std::map<Segment_2, int, Segment2Comparator>
> planar_map(
    const Eigen::MatrixXd &V_2d,
    const Eigen::MatrixXi &E_connectivity
) {
    // Create CGAL Arrangement_2 object
    Arrangement_with_history_2 arr;

    // Store original segments for later reference
    std::vector<std::vector<Point_2>> original_segments;

    // First pass: collect all unique points and store original segments
    std::vector<Segment_2> segments_to_insert;

    std::map<Point_2, int> point_to_node;

    std::map<Segment_2, int, Segment2Comparator> segment_to_segment;

    const auto& N = V_2d;
    const auto& S = E_connectivity;

    std::unordered_map<int, Point_2> point_map;

    for (int j = 0; j < N.rows(); ++j) {
        auto p = Point_2(N(j, 0), N(j, 1));
        point_map.insert({j, p});
        point_to_node.insert({p, j});
    }

    for (int j = 0; j < S.rows(); ++j) {
        int start_idx = S(j, 0);
        int end_idx = S(j, 1);

        Point_2 p1 = point_map[start_idx];
        Point_2 p2 = point_map[end_idx];

        Segment_2 segment(p1, p2);
        segments_to_insert.push_back(segment);
        segment_to_segment.insert({segment, j});
    }

    // Insert all segments into the arrangement
    CGAL::insert(arr, segments_to_insert.begin(), segments_to_insert.end());

    return std::make_tuple(arr, point_to_node, segment_to_segment);
}
}
