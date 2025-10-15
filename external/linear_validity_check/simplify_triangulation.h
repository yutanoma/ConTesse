#pragma once

#include <tuple>
#include <Eigen/Core>

namespace utils {

// Simplify triangulation into polygons; return triangle mesh and mapping
std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXi>
simplify_triangulation(const Eigen::MatrixXd &V,
                       const Eigen::MatrixXi &F);

}
