#include "fast_validity_check.h"

#include "planar_map.h"
#include "assign_qi.h"
#include "validity_check.h"
#include "triangulate_valid_contour.h"

namespace utils {
std::map<Segment_2, int, Segment2Comparator> segment_to_orientation(
    const Arrangement_with_history_2 &arr,
    const std::map<Segment_2, int, Segment2Comparator>
        &segment_to_segment,
    const std::map<Point_2, int> &point_to_node,
    const Eigen::MatrixXd &V_2d,
    const Eigen::MatrixXi &E_connectivity,
    const Eigen::VectorXi &E_sign
) {
  std::map<Segment_2, int, Segment2Comparator> segment_to_orientation;

  for (auto it = arr.edges_begin(); it != arr.edges_end(); it++) {
    Segment_2 this_segment = it->curve();
    Segment_2 original_segment;
    
    if (arr.number_of_originating_curves(it) == 0) {
      original_segment = this_segment;
    } else if (arr.number_of_originating_curves(it) == 1) {
      original_segment = *arr.originating_curves_begin(it);
    } else {
      std::cerr << "Error: it has " << arr.number_of_originating_curves(it)
                << " originating curves" << std::endl;
      throw std::runtime_error("Error: it has " + std::to_string(arr.number_of_originating_curves(it)) + " originating curves");
    }

    if (segment_to_segment.find(original_segment) == segment_to_segment.end()) {
      std::cerr << "Error: original segment not found" << std::endl;
      throw std::runtime_error("Error: original segment not found");
    }

    // get the edge in the original segment
    auto edge_id = segment_to_segment.at(original_segment);
    auto v0 = E_connectivity(edge_id, 0);
    auto v1 = E_connectivity(edge_id, 1);

    Eigen::Vector2d v0_2d(V_2d(v0, 0), V_2d(v0, 1));
    Eigen::Vector2d v1_2d(V_2d(v1, 0), V_2d(v1, 1));

    // get the edge in the new segment
    auto v0_this = this_segment.source();
    auto v1_this = this_segment.target();

    Eigen::Vector2d v0_this_2d(v0_this.x(), v0_this.y());
    Eigen::Vector2d v1_this_2d(v1_this.x(), v1_this.y());

    double d00 = (v0_this_2d - v0_2d).norm();
    double d10 = (v1_this_2d - v0_2d).norm();

    if (d00 < d10) {
      // v0_this is closer to v0_2d than v1_this
      // this means that the orientation matches the original edge
      segment_to_orientation[this_segment] = E_sign(edge_id);
    } else {
      // v1_this is closer to v0_2d than v0_this
      // this means that the orientation is the opposite of the original edge
      segment_to_orientation[this_segment] = -E_sign(edge_id);
    }
  }

  return segment_to_orientation;
}

std::map<Segment_2, bool, Segment2Comparator> segment_to_cut(
  const Arrangement_with_history_2 &arr,
  const std::map<Segment_2, int, Segment2Comparator>
      &segment_to_segment,
  const Eigen::VectorXi &E_is_cut
) {
  std::map<Segment_2, bool, Segment2Comparator> segment_is_cut;

  for (auto it = arr.edges_begin(); it != arr.edges_end(); it++) {
    Segment_2 this_segment = it->curve();
    Segment_2 original_segment;
    
    if (arr.number_of_originating_curves(it) == 0) {
      original_segment = this_segment;
    } else if (arr.number_of_originating_curves(it) == 1) {
      original_segment = *arr.originating_curves_begin(it);
    } else {
      std::cerr << "Error: it has " << arr.number_of_originating_curves(it)
                << " originating curves" << std::endl;
      throw std::runtime_error("Error: it has " + std::to_string(arr.number_of_originating_curves(it)) + " originating curves");
    }

    if (segment_to_segment.find(original_segment) == segment_to_segment.end()) {
      std::cerr << "Error: original segment not found" << std::endl;
      throw std::runtime_error("Error: original segment not found");
    }

    auto edge_id = segment_to_segment.at(original_segment);

    segment_is_cut[this_segment] = E_is_cut(edge_id);
  }
}

std::map<Segment_2, bool, Segment2Comparator>
segment_to_convex(const Arrangement_with_history_2 &arr,
                  const std::map<Point_2, int> &point_to_node,
                  const Eigen::VectorXi &V_is_cusp
) {
  // for each connected component, the segment that has the lowest winding
  // number on the right side face is convex
  // the convex/concave switches on the cusp vertices
}

std::map<Point_2, bool> point_to_singularity(
  const Arrangement_with_history_2 &arr,
  const std::map<Point_2, int> &point_to_node,
  const Eigen::VectorXi &V_is_singularity
) {
  std::map<Point_2, bool> point_is_singularity;

  for (auto it = arr.vertices_begin(); it != arr.vertices_end(); it++) {
    auto point = it->point();

    if (point_to_node.find(point) == point_to_node.end()) {
      point_is_singularity[point] = false;
    } else {
      auto node_id = point_to_node.at(point);
      point_is_singularity[point] = V_is_singularity(node_id);
    }
  }

  return point_is_singularity;
}

std::tuple<std::map<Point_2, Segment_2>, std::map<Point_2, Segment_2>>
segment_to_casing_edges(
    const Arrangement_with_history_2 &arr,
    const std::map<Point_2, int> &point_to_node,
    const Eigen::MatrixXd &V_3d, const Eigen::VectorXi &V_is_cusp,
    const Eigen::Vector3d &camera_pos
) {

}

std::tuple<bool, Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXi>
fast_validity_check(const Eigen::MatrixXd &V_3d, const Eigen::MatrixXd &V_2d,
                    const Eigen::MatrixXi &E_connectivity,
                    const Eigen::VectorXi &E_sign,
                    const Eigen::VectorXi &E_is_cut,
                    const Eigen::VectorXi &V_is_cusp,
                    const Eigen::VectorXi &V_is_singularity,
                    const Eigen::Vector3d &camera_pos) {
  // 1. create the planar map
  // point_to_node is the map from the *original* vertices to the new vertices
  // i.e., does not include any intersection points
  // segment_to_segment is the map from the *original* segments to the new segments
  // i.e., does not include any new segments made by intersections
  auto [arr, point_to_node, segment_to_segment] = planar_map(V_2d, E_connectivity);

  // 2. get the labeling
  auto segment_orientation = segment_to_orientation(arr, segment_to_segment, point_to_node, V_2d, E_connectivity, E_sign);
  auto segment_is_cut = segment_to_cut(arr, segment_to_segment, E_is_cut);
  auto point_is_singularity = point_to_singularity(arr, point_to_node, V_is_singularity);
  auto [upper_casing_edges, lower_casing_edges] = segment_to_casing_edges(arr, point_to_node, V_3d, V_is_cusp, camera_pos);
  auto segment_is_convex = segment_to_convex(arr, point_to_node, V_is_cusp);

  // 3. check the validity
  auto [is_valid, qi, qi_mismatch_positions] = validity_check(
      arr, point_is_singularity, segment_is_convex, segment_is_cut, upper_casing_edges,
      lower_casing_edges, segment_orientation);

  if (!is_valid) {
    return {false, Eigen::MatrixXd(), Eigen::MatrixXi(), Eigen::VectorXi()};
  }

  // 4. tessellate the valid contour
  auto [V_out, F_out, point_to_idx] = triangulate_valid_contour(arr, qi);

  // 5. get the map from the original vertices to the new vertices
  Eigen::VectorXi V_IJ(V_2d.rows());
  {
    for (auto point : point_to_node) {
      auto original_vid = point.second;
      auto new_vid = point_to_idx[point.first];
      V_IJ(original_vid) = new_vid;
    }
  }

  return {true, V_out, F_out, V_IJ};
}
} // namespace utils
