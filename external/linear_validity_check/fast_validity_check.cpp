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

  return segment_is_cut;
}

std::map<Segment_2, bool, Segment2Comparator>
segment_to_convex(const Arrangement_with_history_2 &arr,
                  const std::map<Segment_2, int, Segment2Comparator> &right_wn,
                  const std::map<Point_2, int> &point_to_node,
                  const std::map<Segment_2, bool, Segment2Comparator> &is_cut,
                  const Eigen::VectorXi &V_is_cusp
) {
  // for each connected component, the segment that has the lowest winding
  // number on the right side face is convex
  // the convex/concave switches on the cusp vertices
  std::map<Point_2, bool> point_is_cusp;
  for (auto it = arr.vertices_begin(); it != arr.vertices_end(); it++) {
    auto point = it->point();
    if (point_to_node.find(point) == point_to_node.end()) {
      point_is_cusp[point] = false;
    } else {
      auto node_id = point_to_node.at(point);
      point_is_cusp[point] = V_is_cusp(node_id);
    }
  }

  // first, group all the segments with its connectivity
  std::vector<std::vector<Segment_2>> segments_group;
  std::map<Segment_2, int, Segment2Comparator> segment_to_group;
  std::queue<Halfedge_const_handle> q;
  for (auto it = arr.edges_begin(); it != arr.edges_end(); it++) {
    segment_to_group[it->curve()] = -1;
    q.push(it);
  }

  int added_segments = 0, total_segments = segment_to_group.size();
  int group_idx = 0;
  while (added_segments < total_segments) {
    auto itr = q.front();
    q.pop();

    if (segment_to_group.at(itr->curve()) != -1) {
      continue;
    }

    auto current_itr = itr;
    do {
      auto segment = current_itr->curve();

      if (segment_to_group.at(segment) == -1) {
        segment_to_group[segment] = group_idx;
        segments_group[group_idx].push_back(segment);
        added_segments++;
      } else {
        break;
      }

      // determine the next itr
      auto next_pt = current_itr->target();

      if (next_pt->degree() == 2) {
        current_itr = current_itr->next();
      } else if (next_pt->degree() == 3) {
        // choose the halfedge that is (1) not this itr and (2) not a cut
        auto halfedge = next_pt->incident_halfedges();
        do {
          if (is_cut.at(halfedge->curve())) {
            halfedge++;
            continue;
          }
          if (halfedge == current_itr) {
            halfedge++;
            continue;
          }

          current_itr = halfedge;
          break;
        } while (halfedge != next_pt->incident_halfedges());
      } else if (next_pt->degree() == 4) {
        // choose the halfedge that has the same source segment as the current itr
        auto halfedge = next_pt->incident_halfedges();
        do {
          if (is_cut.at(halfedge->curve())) {
            halfedge++;
            continue;
          }
          if (halfedge == current_itr) {
            halfedge++;
            continue;
          }
          if (arr.number_of_originating_curves(halfedge) != 1) {
            halfedge++;
            continue;
          }
          auto original_segment = *arr.originating_curves_begin(halfedge);
          auto original_segment_source =
              *arr.originating_curves_begin(current_itr);
          if (is_identical_segment(original_segment, original_segment_source)) {
            current_itr = halfedge;
            break;
          } else {
            halfedge++;
            continue;
          }
        } while (halfedge != next_pt->incident_halfedges());
      }
    } while (current_itr != itr);

    group_idx++;
  }

  // second, for each connected component, find the segment that has the lowest
  // winding number on the right side face
  std::map<Segment_2, bool, Segment2Comparator> segment_is_convex;
  std::map<Segment_2, bool, Segment2Comparator> segment_labeled;
  for (auto it = arr.edges_begin(); it != arr.edges_end(); it++) {
    segment_labeled[it->curve()] = false;
  }

  for (int i = 0; i < segments_group.size(); i++) {
    int min_wn = std::numeric_limits<int>::max();
    Segment_2 min_wn_segment;
    for (int j = 0; j < segments_group[i].size(); j++) {
      auto right_side_wn = right_wn.at(segments_group[i][j]);

      if (right_side_wn < min_wn) {
        min_wn = right_side_wn;
        min_wn_segment = segments_group[i][j];
      }
    }

    segment_labeled[min_wn_segment] = true;
    segment_is_convex[min_wn_segment] = true;
  }

  // finally, label the segment as convex or concave
  for (int i = 0; i < segments_group.size(); i++) {
    for (int j = 0; j < 2 * segments_group[i].size(); j++) {
      auto prev_idx = (j + segments_group[i].size() - 1) % segments_group[i].size();
      auto curr_idx = j % segments_group[i].size();

      if (!segment_labeled[segments_group[i][prev_idx]]) {
        continue;
      }

      segment_labeled[segments_group[i][curr_idx]] = true;

      // get the common point
      Point_2 common_point = segments_group[i][prev_idx].source() ==
                                     segments_group[i][curr_idx].source()
                                 ? segments_group[i][prev_idx].source()
                             : segments_group[i][prev_idx].target() ==
                                     segments_group[i][curr_idx].target()
                                 ? segments_group[i][prev_idx].target()
                             : segments_group[i][prev_idx].source() ==
                                     segments_group[i][curr_idx].target()
                                 ? segments_group[i][prev_idx].source()
                             : segments_group[i][prev_idx].target() ==
                                     segments_group[i][curr_idx].source()
                                 ? segments_group[i][prev_idx].target()
                                 : segments_group[i][curr_idx].source();

      if (point_is_cusp.at(common_point)) {
        segment_is_convex[segments_group[i][curr_idx]] = !segment_is_convex[segments_group[i][prev_idx]];
      } else {
        segment_is_convex[segments_group[i][curr_idx]] = segment_is_convex[segments_group[i][prev_idx]];
      }
    }
  }

  return segment_is_convex;
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
  std::map<Point_2, Segment_2> upper_casing_edges;
  std::map<Point_2, Segment_2> lower_casing_edges;

  for (auto it = arr.vertices_begin(); it != arr.vertices_end(); it++) {
    if (it->degree() != 4) {
      continue;
    }

    auto point = it->point();

    auto heit = it->incident_halfedges();
    std::map<Segment_2, std::vector<Segment_2>, Segment2Comparator> original_to_new_segments;
    do {
      Segment_2 this_segment = heit->curve();
      Segment_2 original_segment;

      if (arr.number_of_originating_curves(heit) == 0) {
        original_segment = this_segment;
      } else if (arr.number_of_originating_curves(heit) == 1) {
        original_segment = *arr.originating_curves_begin(heit);
      } else {
        std::cerr << "Error: it has " << arr.number_of_originating_curves(heit) << " originating curves" << std::endl;
        throw std::runtime_error("Error: it has " + std::to_string(arr.number_of_originating_curves(heit)) + " originating curves");
      }

      original_to_new_segments[original_segment].push_back(this_segment);

      ++heit;
    } while (heit != it->incident_halfedges());

    if (original_to_new_segments.size() != 2) {
      std::cerr << "Error: original_to_new_segments.size() != 2" << std::endl;
      throw std::runtime_error("Error: original_to_new_segments.size() != 2");
    }

    Segment_2 original_segment_1 = original_to_new_segments.begin()->first;
    Segment_2 original_segment_2 = original_to_new_segments.begin()++->first;

    Point_2 original_point_1_1 = original_segment_1.source();
    Point_2 original_point_1_2 = original_segment_1.target();
    Point_2 original_point_2_1 = original_segment_2.source();
    Point_2 original_point_2_2 = original_segment_2.target();

    Eigen::Vector2d original_point_1_1_2d(original_point_1_1.x(), original_point_1_1.y());
    Eigen::Vector2d original_point_1_2_2d(original_point_1_2.x(), original_point_1_2.y());
    Eigen::Vector2d original_point_2_1_2d(original_point_2_1.x(), original_point_2_1.y());
    Eigen::Vector2d original_point_2_2_2d(original_point_2_2.x(),
                                          original_point_2_2.y());

    Eigen::Vector2d intersection_2d(
        point.x(), point.y()
    );

    double d_1_total = (original_point_1_1_2d - original_point_1_2_2d).norm();
    double d_2_total = (original_point_2_1_2d - original_point_2_2_2d).norm();

    double d_1_int = (original_point_1_1_2d - intersection_2d).norm();
    double d_2_int = (original_point_2_1_2d - intersection_2d).norm();

    int p1_1_idx = point_to_node.at(original_point_1_1);
    int p1_2_idx = point_to_node.at(original_point_1_2);
    int p2_1_idx = point_to_node.at(original_point_2_1);
    int p2_2_idx = point_to_node.at(original_point_2_2);

    Eigen::Vector3d p1_1_3d(V_3d(p1_1_idx, 0), V_3d(p1_1_idx, 1), V_3d(p1_1_idx, 2));
    Eigen::Vector3d p1_2_3d(V_3d(p1_2_idx, 0), V_3d(p1_2_idx, 1), V_3d(p1_2_idx, 2));
    Eigen::Vector3d p2_1_3d(V_3d(p2_1_idx, 0), V_3d(p2_1_idx, 1), V_3d(p2_1_idx, 2));
    Eigen::Vector3d p2_2_3d(V_3d(p2_2_idx, 0), V_3d(p2_2_idx, 1), V_3d(p2_2_idx, 2));

    Eigen::Vector3d p1_int = p1_1_3d + (p1_2_3d - p1_1_3d) * d_1_int / d_1_total;
    Eigen::Vector3d p2_int = p2_1_3d + (p2_2_3d - p2_1_3d) * d_2_int / d_2_total;

    double p1_int_camera = (p1_int - camera_pos).norm();
    double p2_int_camera = (p2_int - camera_pos).norm();

    Segment_2 upper_segment, lower_segment;

    if (p1_int_camera < p2_int_camera) {
      upper_segment = original_segment_1;
      lower_segment = original_segment_2;
    } else {
      upper_segment = original_segment_2;
      lower_segment = original_segment_1;
    }

    upper_casing_edges[point] = upper_segment;
    lower_casing_edges[point] = lower_segment;
  }

  return {upper_casing_edges, lower_casing_edges};
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
  auto [arr, point_to_node, segment_to_segment] =
      planar_map(V_2d, E_connectivity);

  // 2. get the labeling
  auto segment_orientation = segment_to_orientation(arr, segment_to_segment, point_to_node, V_2d, E_connectivity, E_sign);
  auto segment_is_cut = segment_to_cut(arr, segment_to_segment, E_is_cut);
  auto point_is_singularity = point_to_singularity(arr, point_to_node, V_is_singularity);
  auto [upper_casing_edges, lower_casing_edges] =
      segment_to_casing_edges(arr, point_to_node, V_3d, V_is_cusp, camera_pos);

  auto wn = face_wn(arr, segment_orientation);
  auto [left_face, right_face] = segment_to_face(arr, segment_orientation);
  std::map<Segment_2, int, Segment2Comparator> right_wn;
  for (auto it = arr.edges_begin(); it != arr.edges_end(); it++) {
    auto segment = it->curve();
    auto rf = right_face.at(segment);
    right_wn[segment] = wn.at(rf);
  }

  auto segment_is_convex = segment_to_convex(arr, right_wn, point_to_node, segment_is_cut, V_is_cusp);

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
