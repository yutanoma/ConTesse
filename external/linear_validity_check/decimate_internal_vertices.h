#pragma once

#include <tuple>
#include <Eigen/Core>

namespace utils {

// Remove internal vertices via edge flips while preserving marked vertices
std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXi>
decimate_internal_vertices(const Eigen::MatrixXd &V,
                           const Eigen::MatrixXi &F,
                           const Eigen::VectorXi &is_decimating);

}


