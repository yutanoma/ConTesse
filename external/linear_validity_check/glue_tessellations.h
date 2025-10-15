#pragma once

#include "planar_map.h"

namespace utils {

// tessellation of a face
// includes:
// - V: vertices
// - F: faces
// - EV: edges of the face
// - FE: faces of the edges
// - EF: edges of the faces
// - point_to_vid: map from point to vertex id
// - segment_to_eid: map from segment to edge id
using tessellation_t =
    std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXi,
               Eigen::MatrixXi, Eigen::MatrixXi, std::map<Point_2, int>,
               std::map<Segment_2, int, Segment2Comparator>>;

// glue tessellations across faces/layers into a single vertex/face list
std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, std::vector<Point_2>>
glue_tessellations(
    const std::map<std::pair<int, Face_const_handle>, tessellation_t> &face_to_tessellation,
    const std::map<std::pair<int, Segment_2>, bool, IntSegmentPairComparator> &left_side_exists,
    const std::map<std::pair<int, Segment_2>, bool, IntSegmentPairComparator> &right_side_exists,
    const std::map<std::pair<int, Segment_2>, std::pair<int, Face_const_handle>, IntSegmentPairComparator> &left_side_face,
    const std::map<std::pair<int, Segment_2>, std::pair<int, Face_const_handle>, IntSegmentPairComparator> &right_side_face);
} // namespace utils
