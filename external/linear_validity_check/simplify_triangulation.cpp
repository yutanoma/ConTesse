#include "simplify_triangulation.h"

namespace utils {

std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXi>
simplify_triangulation(const Eigen::MatrixXd &V,
                       const Eigen::MatrixXi &F)
{
  Eigen::VectorXi V_simp_to_out;
  V_simp_to_out.setLinSpaced(V.rows(), 0, V.rows() - 1);
  return {V, F, V_simp_to_out};
}
}
