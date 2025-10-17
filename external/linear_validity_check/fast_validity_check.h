#pragma once

#include "planar_map.h"

namespace utils {
std::tuple<bool, Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXi>
fast_validity_check(const Eigen::MatrixXd &V_3d, const Eigen::MatrixXd &V_2d,
                    const Eigen::MatrixXi &E_connectivity,
                    const Eigen::VectorXi &E_sign,
                    const Eigen::VectorXi &E_is_cut,
                    const Eigen::VectorXi &V_is_cusp,
                    const Eigen::VectorXi &V_is_singularity,
                    const Eigen::Vector3d &camera_pos,
                    bool is_back_facing,
                    bool is_debug = false);
} // namespace utils