#pragma once

#include "planar_map.h"

namespace utils {
// 1. triangulate each face w/ winding number > 0 in the arrangement
// 2. connect the vertices with the connected components
// 3. iteratively remove the internal vertices via edge flips (decimate_internal_vertices.h)
// 4. divide it into simple polygons using Shor van Wyck's `simplify_triangulation` function (simplify_polygon.h)
std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, std::vector<Point_2>>
triangulate_valid_contour(
    const Arrangement_with_history_2 &arr,
    const std::map<Segment_2, int, Segment2Comparator> &qi,
    const std::map<Segment_2, int, Segment2Comparator> &segment_orientation);
} // namespace utils
