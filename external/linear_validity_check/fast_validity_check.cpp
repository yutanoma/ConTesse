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

    Eigen::Vector2d v0_this_2d(CGAL::to_double(v0_this.x()), CGAL::to_double(v0_this.y()));
    Eigen::Vector2d v1_this_2d(CGAL::to_double(v1_this.x()), CGAL::to_double(v1_this.y()));

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
                  const std::map<Segment_2, int, Segment2Comparator> &segment_to_segment,
                  const Eigen::MatrixXi &E_connectivity,
                  const Eigen::VectorXi &E_is_cut,
                  const Eigen::VectorXi &V_is_cusp
) {
  // for each connected component, the segment that has the lowest winding
  // number on the right side face is convex
  // the convex/concave switches on the cusp vertices

  // 1. for each segment, find the minimum winding number on the right side face
  Eigen::VectorXi right_min_wn = Eigen::VectorXi::Constant(E_connectivity.rows(), std::numeric_limits<int>::max());
  for (auto it = arr.edges_begin(); it != arr.edges_end(); it++) {
    auto segment = it->curve();
    auto right_face = right_wn.at(segment);

    Segment_2 original_segment;
    if (arr.number_of_originating_curves(it) == 0) {
      original_segment = segment;
    } else if (arr.number_of_originating_curves(it) == 1) {
      original_segment = *arr.originating_curves_begin(it);
    } else {
      std::cerr << "Error: segment has " << arr.number_of_originating_curves(it) << " originating curves" << std::endl;
      throw std::runtime_error("Error: segment has " + std::to_string(arr.number_of_originating_curves(it)) + " originating curves");
    }

    auto edge_id = segment_to_segment.at(original_segment);

    if (right_min_wn(edge_id) > right_face) {
      right_min_wn(edge_id) = right_face;
    }
  }

  // 2. divide the segments into patches
  std::vector<std::vector<int>> vertex_to_edge(V_is_cusp.rows());
  for (int i = 0; i < E_connectivity.rows(); i++) {
    if (!E_is_cut(i)) {
      vertex_to_edge[E_connectivity(i, 0)].push_back(i);
      vertex_to_edge[E_connectivity(i, 1)].push_back(i);
    }
  }

  for (int i = 0; i < vertex_to_edge.size(); i++) {
    if (vertex_to_edge[i].size() !=  0 && vertex_to_edge[i].size() != 2) {
      std::cerr << "Error: vertex_to_edge[" << i << "] has "
                << vertex_to_edge[i].size() << " edges" << std::endl;
      std::cout << "This means that the loop is not a manifold" << std::endl;
      throw std::runtime_error("Error: vertex_to_edge[" + std::to_string(i) + "] has " + std::to_string(vertex_to_edge[i].size()) + " edges");
    }
  }

  std::vector<std::vector<int>> edge_groups;
  Eigen::VectorXi edge_to_group = Eigen::VectorXi::Constant(E_connectivity.rows(), -1);
  for (int i = 0; i < E_connectivity.rows(); i++) {
    if (E_is_cut(i)) {
      continue;
    }

    if (edge_to_group(i) != -1) {
      continue;
    }

    if (vertex_to_edge[E_connectivity(i, 0)].size() < 2) {
      continue;
    }

    int group_id = edge_groups.size();
    std::vector<int> group;
    int current_edge_id = i;

    while (true) {
      group.push_back(current_edge_id);
      edge_to_group(current_edge_id) = group_id;

      std::vector<int> next_edges;
      for (int j = 0; j < 2; j++) {
        auto vidx = E_connectivity(current_edge_id, j);
        for (auto edge_id : vertex_to_edge[vidx]) {
          if (edge_to_group(edge_id) == -1 && edge_id != current_edge_id) {
            next_edges.push_back(edge_id);
          }
        }
      }

      if (next_edges.empty()) {
        break;
      }

      current_edge_id = next_edges[0];
    }

    edge_groups.push_back(group);
  }

  // 3. for each group, find the minimum winding number
  std::vector<int> group_min_wn_edge(edge_groups.size(), -1);
  for (int i = 0; i < edge_groups.size(); i++) {
    int min_wn = std::numeric_limits<int>::max();
    for (auto edge_id : edge_groups[i]) {
      if (right_min_wn(edge_id) < min_wn) {
        min_wn = right_min_wn(edge_id);
        group_min_wn_edge[i] = edge_id;
      }
    }
  }

  // 4. for each group, flip the convex/concave labeling
  Eigen::VectorXi segment_is_convex = Eigen::VectorXi::Constant(E_connectivity.rows(), -1);
  for (int i = 0; i < edge_groups.size(); i++) {
    int edge_id = group_min_wn_edge[i];
    segment_is_convex(edge_id) = 1;
  }

  for (int i = 0; i < edge_groups.size(); i++) {
    for (int j = 0; j < 2 * edge_groups[i].size(); j++) {
      int edge_id = edge_groups[i][j % edge_groups[i].size()];
      int next_edge_id = edge_groups[i][(j + 1) % edge_groups[i].size()];

      int shared_vertex = -1;
      for (int k = 0; k < 2 && shared_vertex == -1; k++) {
        for (int l = 0; l < 2; l++) {
          if (E_connectivity(edge_id, k) == E_connectivity(next_edge_id, l)) {
            shared_vertex = E_connectivity(edge_id, k);
            break;
          }
        }
      }

      if (shared_vertex == -1) {
        std::cerr << "Error: shared vertex not found" << std::endl;
        std::cout << "edge_id: " << edge_id << std::endl;
        std::cout << "next_edge_id: " << next_edge_id << std::endl;
        std::cout << "E_connectivity(edge_id, 0): " << E_connectivity(edge_id, 0) << std::endl;
        std::cout << "E_connectivity(edge_id, 1): " << E_connectivity(edge_id, 1) << std::endl;
        std::cout << "E_connectivity(next_edge_id, 0): " << E_connectivity(next_edge_id, 0) << std::endl;
        std::cout << "E_connectivity(next_edge_id, 1): " << E_connectivity(next_edge_id, 1) << std::endl;
        throw std::runtime_error("Error: shared vertex not found");
      }

      if (segment_is_convex(edge_id) == -1) {
        continue;
      }

      int prev_label = segment_is_convex(edge_id);
      int next_label = V_is_cusp(shared_vertex) ? !prev_label : prev_label;

      segment_is_convex(next_edge_id) = next_label;
    }
  }

  // 5. return the result
  std::map<Segment_2, bool, Segment2Comparator> segment_is_convex_result;
  for (auto it = arr.edges_begin(); it != arr.edges_end(); it++) {
    Segment_2 segment = it->curve();
    Segment_2 original_segment;

    if (arr.number_of_originating_curves(it) == 0) {
      original_segment = segment;
    } else if (arr.number_of_originating_curves(it) == 1) {
      original_segment = *arr.originating_curves_begin(it);
    } else {
      std::cerr << "Error: segment has " << arr.number_of_originating_curves(it) << " originating curves" << std::endl;
      throw std::runtime_error("Error: segment has " + std::to_string(arr.number_of_originating_curves(it)) + " originating curves");
    }

    segment_is_convex_result[segment] = segment_is_convex(segment_to_segment.at(original_segment));
  }

  return segment_is_convex_result;
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

std::tuple<std::map<Point_2, std::vector<Segment_2>>, std::map<Point_2, std::vector<Segment_2>>>
segment_to_casing_edges(
    const Arrangement_with_history_2 &arr,
    const std::map<Point_2, int> &point_to_node,
    const Eigen::MatrixXd &V_3d, const Eigen::VectorXi &V_is_cusp,
    const Eigen::Vector3d &camera_pos
) {
  std::map<Point_2, std::vector<Segment_2>> upper_casing_edges;
  std::map<Point_2, std::vector<Segment_2>> lower_casing_edges;

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

    // handle degenerate cases
    if (original_to_new_segments.size() != 2) {
      // group it with the segments that are connected
      std::map<Segment_2, bool, Segment2Comparator> segment_is_grouped;
      std::map<Segment_2, std::vector<Segment_2>, Segment2Comparator> _original_to_new_segments;
      for (auto it0 : original_to_new_segments) {
        for (auto it1 : original_to_new_segments) {
          if (segment_is_grouped.find(it0.first) != segment_is_grouped.end() ||
              segment_is_grouped.find(it1.first) != segment_is_grouped.end()
          ) {
            continue;
          }

          if (is_identical_segment(it0.first, it1.first)) {
            continue;
          }

          auto segment_0 = it0.first;
          auto segment_1 = it1.first;

          if (segment_0.source() == segment_1.source() ||
              segment_0.source() == segment_1.target() ||
              segment_0.target() == segment_1.source() ||
              segment_0.target() == segment_1.target()) {
            auto segments_0 = it0.second;
            auto segments_1 = it1.second;

            segment_is_grouped[segment_0] = true;
            segment_is_grouped[segment_1] = true;

            _original_to_new_segments[segment_0].insert(
                _original_to_new_segments[segment_0].end(), segments_0.begin(),
                segments_0.end());
            _original_to_new_segments[segment_0].insert(
                _original_to_new_segments[segment_0].end(), segments_1.begin(),
                segments_1.end());

            break;
          }
        }

        if (segment_is_grouped.find(it0.first) == segment_is_grouped.end()) {
          _original_to_new_segments[it0.first] = it0.second;
          segment_is_grouped[it0.first] = true;
        }
      }

      original_to_new_segments = _original_to_new_segments;
    }

    if (
      original_to_new_segments.size() != 2 ||
      original_to_new_segments.begin()->second.size() != 2 ||
      (original_to_new_segments.begin()++)->second.size() != 2
    ) {
      std::cerr << "Error: original_to_new_segments.size() != 2" << std::endl;
      std::cerr << "original_to_new_segments.size(): " << original_to_new_segments.size() << std::endl;
      for (auto it : original_to_new_segments) {
        std::cout << "original_segment: " << it.first << std::endl;
        for (auto s : it.second) {
          std::cout << "   this_segment: " << s << std::endl;
        }
        std::cout << std::endl;
      }
      throw std::runtime_error("Error: original_to_new_segments.size() != 2");
    }

    Segment_2 original_segment_1 = original_to_new_segments.begin()->first;
    Segment_2 original_segment_2 = original_to_new_segments.begin()++->first;

    Point_2 original_point_1_1 = original_segment_1.source();
    Point_2 original_point_1_2 = original_segment_1.target();
    Point_2 original_point_2_1 = original_segment_2.source();
    Point_2 original_point_2_2 = original_segment_2.target();

    Eigen::Vector2d original_point_1_1_2d(CGAL::to_double(original_point_1_1.x()), CGAL::to_double(original_point_1_1.y()));
    Eigen::Vector2d original_point_1_2_2d(CGAL::to_double(original_point_1_2.x()), CGAL::to_double(original_point_1_2.y()));
    Eigen::Vector2d original_point_2_1_2d(CGAL::to_double(original_point_2_1.x()), CGAL::to_double(original_point_2_1.y()));
    Eigen::Vector2d original_point_2_2_2d(CGAL::to_double(original_point_2_2.x()),
                                          CGAL::to_double(original_point_2_2.y()));

    Eigen::Vector2d intersection_2d(
        CGAL::to_double(point.x()), CGAL::to_double(point.y())
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

    upper_casing_edges[point] = original_to_new_segments[upper_segment];
    lower_casing_edges[point] = original_to_new_segments[lower_segment];
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
  // segment_to_segment is the map from the *original* segments to the new
  // segments i.e., does not include any new segments made by intersections
  std::cout << "constructing planar map" << std::endl;
  auto [arr, point_to_node, segment_to_segment] =
      planar_map(V_2d, E_connectivity);
  std::cout << "constructed planar map" << std::endl;

  // 2. get the labeling
  std::cout << "getting labeling" << std::endl;
  auto segment_orientation = segment_to_orientation(
      arr, segment_to_segment, point_to_node, V_2d, E_connectivity, E_sign);
  std::cout << "got segment orientation" << std::endl;

  std::cout << "getting segment is cut" << std::endl;
  auto segment_is_cut = segment_to_cut(arr, segment_to_segment, E_is_cut);
  std::cout << "got segment is cut" << std::endl;

  std::cout << "getting point is singularity" << std::endl;
  auto point_is_singularity = point_to_singularity(arr, point_to_node, V_is_singularity);
  std::cout << "got point is singularity" << std::endl;

  std::cout << "getting upper and lower casing edges" << std::endl;
  auto [upper_casing_edges, lower_casing_edges] =
      segment_to_casing_edges(arr, point_to_node, V_3d, V_is_cusp, camera_pos);
  std::cout << "got upper and lower casing edges" << std::endl;

  std::cout << "getting face winding numbers" << std::endl;
  auto wn = face_wn(arr, segment_orientation);
  std::cout << "got face winding numbers" << std::endl;

  std::cout << "getting left and right faces" << std::endl;
  auto [left_face, right_face] = segment_to_face(arr, segment_orientation);
  std::map<Segment_2, int, Segment2Comparator> right_wn;
  for (auto it = arr.edges_begin(); it != arr.edges_end(); it++) {
    auto segment = it->curve();
    auto rf = right_face.at(segment);
    right_wn[segment] = wn.at(rf);
  }
  std::cout << "got left and right faces" << std::endl;

  std::cout << "getting segment is convex" << std::endl;
  auto segment_is_convex = segment_to_convex(arr, right_wn, point_to_node,
                                             segment_to_segment, E_connectivity,
                                             E_is_cut, V_is_cusp);
  std::cout << "got segment is convex" << std::endl;

  // 3. check the validity
  std::cout << "checking validity" << std::endl;
  auto [is_valid, qi, qi_mismatch_positions] = validity_check(
      arr, point_is_singularity, segment_is_convex, segment_is_cut, upper_casing_edges,
      lower_casing_edges, segment_orientation);
  std::cout << "checked validity" << std::endl;

  if (!is_valid) {
    return {false, Eigen::MatrixXd(), Eigen::MatrixXi(), Eigen::VectorXi()};
  }

  // 4. tessellate the valid contour
  std::cout << "tessellating valid contour" << std::endl;
  auto [V_out, F_out, V_to_point] = triangulate_valid_contour(arr, qi, segment_orientation, segment_is_cut);
  std::cout << "tessellated valid contour" << std::endl;

  // 5. get the map from the original vertices to the new vertices
  Eigen::VectorXi V_JI(V_out.rows());
  {
    for (int vid_final = 0; vid_final < V_to_point.size(); vid_final++) {
      auto point = V_to_point[vid_final];

      if (point_to_node.find(point) == point_to_node.end()) {
        std::cerr << "Error: point not found" << std::endl;
        throw std::runtime_error("Error: point not found");
      }

      auto original_vid = point_to_node.at(point);
      V_JI(vid_final) = original_vid;
    }
  }

  return {true, V_out, F_out, V_JI};
}
} // namespace utils
