#pragma once

#include "planar_map.h"

namespace utils {

// Triangulate a face and return tessellation data:
// V (|V|x2), F (|F|x3), EV, FE, EF adjacency, and point/segment maps
std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXi,
           Eigen::MatrixXi, Eigen::MatrixXi, std::map<Point_2, int>,
           std::map<Segment_2, int, Segment2Comparator>>
tessellate_face(const Arrangement_with_history_2 &arr,
                const Face_const_handle &face);

}
