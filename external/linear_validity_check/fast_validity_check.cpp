#include "fast_validity_check.h"

namespace utils {
std::tuple<bool, Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXi>
fast_validity_check(const Eigen::MatrixXd &V_3d, const Eigen::MatrixXd &V_2d,
                    const Eigen::MatrixXi &E_connectivity,
                    const Eigen::VectorXi &E_sign,
                    const Eigen::VectorXi &E_is_cut,
                    const Eigen::VectorXi &V_is_cusp,
                    const Eigen::VectorXi &V_is_singularity,
                    const Eigen::Vector3d &camera_pos) {
  // 1. create the planar map
  std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>> contours_2d = {
    {V_2d, E_connectivity}
  };

  auto [arr, point_to_node, segment_to_segment] = planar_map(contours_2d);

  // 2. get the labeling
  auto segment_orientation = segment_to_orientation(arr, segment_to_segment, E_connectivity, E_sign);
  auto segment_is_convex = segment_to_convex(arr, point_to_node, V_is_cusp);
  auto segment_is_cut = segment_to_cut(arr, segment_to_segment, E_is_cut);
  auto [upper_casing_edges, lower_casing_edges] = segment_to_casing_edges(arr, point_to_node, V_3d, V_is_cusp, camera_pos);

  // 3. check the validity
  auto [is_valid, qi, qi_mismatch_positions] = validity_check(
      arr, point_to_node, segment_is_convex, segment_is_cut, upper_casing_edges,
      lower_casing_edges, segment_orientation);

  if (!is_valid) {
    return {false, Eigen::MatrixXd(), Eigen::MatrixXi(), Eigen::VectorXi()};
  }

  // 4. tessellate the valid contour
  auto [V_out, F_out] = tessellate_valid_contour(
      arr, point_to_node, segment_is_convex, segment_is_cut, upper_casing_edges,
      lower_casing_edges, segment_orientation, qi);

  return {true, V_out, F_out, qi_mismatch_positions};
}
} // namespace utils
