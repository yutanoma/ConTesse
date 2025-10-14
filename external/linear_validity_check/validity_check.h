#pragma once

#include "planar_map.h"

namespace utils {
std::tuple<
    std::map<Segment_2, Face_const_handle, Segment2Comparator>, // left face
    std::map<Segment_2, Face_const_handle, Segment2Comparator>  // right face
> segment_to_face(const Arrangement_with_history_2 &arr,
  const std::map<Segment_2, int, Segment2Comparator> &segment_orientation);

std::map<Face_const_handle, int> face_wn(const Arrangement_with_history_2 &arr);

// check if the qi is consistent with the arrangement
// it checks if
// (1) the qi matches the topological constraints of an image-space intersection, and
// (2) the qi matches the winding number constraints
// no need for the convexity or casing information
std::tuple<bool, std::vector<Point_2>>
check_qi_mismatch(const Arrangement_with_history_2 &arr,
                  const std::map<Segment_2, int, Segment2Comparator> &segment_orientation,
                  const std::map<Segment_2, int, Segment2Comparator> &qi,
                  // invalid segments will be ignored for this check
                  const std::map<Segment_2, bool, Segment2Comparator> &is_segment_invalid
);
} // namespace utils
