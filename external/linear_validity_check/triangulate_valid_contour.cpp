#include "triangulate_valid_contour.h"
#include "validity_check.h"
#include "tessellate_face.h"
#include "decimate_internal_vertices.h"
#include "simplify_triangulation.h"
#include "glue_tessellations.h"

namespace utils {
std::tuple <
  std::map<std::pair<int, Segment_2>, bool, IntSegmentPairComparator>,
  std::map<std::pair<int, Segment_2>, bool, IntSegmentPairComparator>,
  std::map<std::pair<int, Segment_2>, std::pair<int, Face_const_handle>, IntSegmentPairComparator>,
  std::map<std::pair<int, Segment_2>, std::pair<int, Face_const_handle>, IntSegmentPairComparator>
> get_left_right_side_faces(
  const Arrangement_with_history_2 &arr,
  const std::map<Face_const_handle, int> &wn,
  const std::map<Segment_2, int, Segment2Comparator> &qi,
  const std::map<Segment_2, int, Segment2Comparator> &segment_orientation,
  const std::map<Segment_2, bool, Segment2Comparator> &segment_is_cut
) {
  std::map<std::pair<int, Segment_2>, bool, IntSegmentPairComparator> left_side_exists;
  std::map<std::pair<int, Segment_2>, bool, IntSegmentPairComparator> right_side_exists;

  std::map<std::pair<int, Segment_2>, std::pair<int, Face_const_handle>, IntSegmentPairComparator>
      right_side_face;
  std::map<std::pair<int, Segment_2>, std::pair<int, Face_const_handle>, IntSegmentPairComparator>
      left_side_face;

  auto [segment_to_left_face, segment_to_right_face] = segment_to_face(arr, segment_orientation);

  for (auto it = arr.edges_begin(); it != arr.edges_end(); it++) {
    auto segment = it->curve();

    // get the wn on the left side and right side
    auto left_face = segment_to_left_face.at(segment);
    auto right_face = segment_to_right_face.at(segment);
    auto left_wn = wn.at(left_face);
    auto right_wn = wn.at(right_face);
    (void)left_wn; // left_wn is currently unused

    if (segment_is_cut.at(segment)) {
      for (int i = 0; i < right_wn + 1; i++) {
        left_side_exists[std::make_pair(i, segment)] = false;
        right_side_exists[std::make_pair(i, segment)] = false;
      }
      continue;
    }

    int layers_num = right_wn + 1;
    int edge_qi = qi.at(segment);

    for (int i = 0; i < layers_num; i++) {
      if (i == edge_qi) {
        left_side_exists[std::make_pair(i, segment)] = true;
        right_side_exists[std::make_pair(i, segment)] = false;
      } else {
        left_side_exists[std::make_pair(i, segment)] = true;
        right_side_exists[std::make_pair(i, segment)] = true;
      }
    }

    for (int i = 0; i < edge_qi; i++) {
      // for each layer above the boundary edge, the layers of the same qi are
      // connected
      left_side_face[std::make_pair(i, segment)] = std::make_pair(i, left_face);
      right_side_face[std::make_pair(i, segment)] = std::make_pair(i, right_face);
    }

    for (int i = edge_qi + 1; i < layers_num; i++) {
      // for each layer below the boundary edge, the layers of the different qi
      // are connected
      left_side_face[std::make_pair(i, segment)] = std::make_pair(i, left_face);
      right_side_face[std::make_pair(i, segment)] = std::make_pair(i - 1, right_face);
    }

    // the layer of the boundary edge
    left_side_face[std::make_pair(edge_qi, segment)] = std::make_pair(edge_qi, left_face);
  }

  return {left_side_exists, right_side_exists, left_side_face, right_side_face};
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, std::vector<Point_2>>
triangulate_valid_contour(
    const Arrangement_with_history_2 &arr,
    const std::map<Segment_2, int, Segment2Comparator> &qi,
    const std::map<Segment_2, int, Segment2Comparator> &segment_orientation,
    const std::map<Segment_2, bool, Segment2Comparator> &segment_is_cut) {
  auto wn = face_wn(arr, segment_orientation);

  // 1. get the left and right side faces
  auto [left_side_exists, right_side_exists, left_side_face, right_side_face] = get_left_right_side_faces(arr, wn, qi, segment_orientation, segment_is_cut);

  // 2. get the tessellation of each faces
  std::map<std::pair<int, Face_const_handle>, tessellation_t> face_to_tessellation;
  for (auto it = arr.faces_begin(); it != arr.faces_end(); it++) {
    if (wn[it] <= 0) {
      continue;
    }

    auto tessellation = tessellate_face(arr, it);
    int wn_face = wn.at(it);

    for (int i = 0; i < wn_face; i++) {
      face_to_tessellation[std::make_pair(i, it)] = tessellation;
    }
  }

  // 3. for each left and right pair, glue the vertices
  auto [V, F, V_to_point] = glue_tessellations(face_to_tessellation, left_side_exists, right_side_exists, left_side_face, right_side_face);

  // 4. remove the internal vertices
  std::map<Point_2, int> point_valences;
  for (auto it = arr.vertices_begin(); it != arr.vertices_end(); it++) {
    auto point = it->point();
    point_valences[point] = it->degree();
  }
  Eigen::VectorXi is_image_intersection(V.rows());
  for (int i = 0; i < V.rows(); i++) {
    auto point = V_to_point[i];
    if (point_valences[point] == 4) {
      is_image_intersection(i) = true;
    } else {
      is_image_intersection(i) = false;
    }
  }
  auto [V_out, F_out, V_out_to_V] = decimate_internal_vertices(V, F, is_image_intersection);

  // 5. simplify the triangulation
  auto [V_simp, F_simp, V_simp_to_out] = simplify_triangulation(V_out, F_out);

  // 6. return the result
  std::vector<Point_2> point_to_idx_final(V_simp.rows());
  for (int i = 0; i < V_simp.rows(); i++) {
    int id_simp = i;
    int id_out = V_simp_to_out(id_simp);
    int id_v = V_out_to_V(id_out);
    Point_2 point = V_to_point[id_v];

    point_to_idx_final[id_simp] = point;
  }

  return {V_simp, F_simp, point_to_idx_final};
}
} // namespace utils
