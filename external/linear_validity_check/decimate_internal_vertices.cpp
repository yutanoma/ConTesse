#include "decimate_internal_vertices.h"

namespace utils {

std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXi>
decimate_internal_vertices(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                           const Eigen::VectorXi &is_decimating
) {
  // remove if the vertex is internal or a specified vertex on the boundary

  // (edge flip the one ring of internal vertices so that it becomes a 3-valence
  // vertex (a 2-valence vertex if on the boundary), which can be trivially removed)


  Eigen::VectorXi V_out_to_V;
  V_out_to_V.setLinSpaced(V.rows(), 0, V.rows() - 1);
  return {V, F, V_out_to_V};
}

}
