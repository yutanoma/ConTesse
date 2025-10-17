#pragma once

#include "planar_map.h"
#include <map>
#include <vector>
#include <string>

namespace utils {
std::string planar_map_to_svg(
    const Arrangement_with_history_2 &arr,
    const std::map<Segment_2, int, Segment2Comparator> &qi,
    const std::map<Segment_2, bool, Segment2Comparator> &segment_is_convex,
    const std::map<Point_2, bool> &point_is_singularity,
    const std::vector<Point_2> &invalid_points
);
}
