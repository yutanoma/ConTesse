#include "glue_tessellations.h"
#include <igl/remove_unreferenced.h>

namespace utils {
// unify the tessellations into a single mesh
// return:
// - face_to_tessellation_unified: map from (qi, face) to the list of face ids in the global tessellation. 
// - V_unified: unified vertices
// - F_unified: unified faces
std::tuple <
    std::map<std::pair<int, Face_const_handle>, Eigen::VectorXi>,
    Eigen::MatrixXd,
    Eigen::MatrixXi,
    std::vector<Point_2>
>
unify_tessellations(
    const std::map<std::pair<int, Face_const_handle>, tessellation_t>
        &face_to_tessellation) {
  int vnum = 0, fnum = 0;
  for (auto it = face_to_tessellation.begin(); it != face_to_tessellation.end();
       it++) {
    auto [V, F, EV, FE, EF, point_to_vid, segment_to_eid] = it->second;
    vnum += V.rows();
    fnum += F.rows();
  }

  Eigen::MatrixXd V_unified(vnum, 2);
  Eigen::MatrixXi F_unified(fnum, 3);
  std::map<std::pair<int, Face_const_handle>, Eigen::VectorXi> face_to_tessellation_unified;
  std::vector<Point_2> V_unified_to_point(vnum);
  int vidx = 0, fidx = 0;
  for (auto it = face_to_tessellation.begin(); it != face_to_tessellation.end();
       it++) {
    auto [V, F, EV, FE, EF, point_to_vid, segment_to_eid] = it->second;
    V_unified.block(vidx, 0, V.rows(), 2) = V;

    F += Eigen::MatrixXi::Constant(F.rows(), 3, vidx);
    F_unified.block(fidx, 0, F.rows(), 3) = F;

    Eigen::VectorXi face_ids(F.rows());
    for (int i = 0; i < F.rows(); i++) {
      face_ids(i) = fidx + i;
    }
    face_to_tessellation_unified[it->first] = face_ids;

    for (auto it : point_to_vid) {
      V_unified_to_point[it.second] = it.first;
    }

    vidx += V.rows();
    fidx += F.rows();
  }

  return {face_to_tessellation_unified, V_unified, F_unified, V_unified_to_point};
}

std::vector<std::pair<int, int>> get_glueing_pairs(
    const std::map<std::pair<int, Face_const_handle>, tessellation_t>
        &face_to_tessellation,
    const std::map<std::pair<int, Face_const_handle>, Eigen::VectorXi>
        &face_to_tessellation_unified,
    const std::map<std::pair<int, Segment_2>, bool, IntSegmentPairComparator>
        &left_side_exists,
    const std::map<std::pair<int, Segment_2>, bool, IntSegmentPairComparator>
        &right_side_exists,
    const std::map<std::pair<int, Segment_2>, std::pair<int, Face_const_handle>,
                   IntSegmentPairComparator> &left_side_face,
    const std::map<std::pair<int, Segment_2>, std::pair<int, Face_const_handle>,
                   IntSegmentPairComparator> &right_side_face) {
  std::map<std::pair<int, Face_const_handle>, Eigen::MatrixXi> face_to_local_EF;
  std::map<std::pair<int, Face_const_handle>, std::map<Segment_2, int, Segment2Comparator>> face_to_local_segment_to_eid;
  for (auto it : face_to_tessellation) {
    auto [V, F, EV, FE, EF, point_to_vid, segment_to_eid] = it.second;
    face_to_local_EF[it.first] = EF;
    face_to_local_segment_to_eid[it.first] = segment_to_eid;
  }

  std::vector<std::pair<int, int>> glueing_pairs;
  for (auto it : left_side_exists) {
    auto lse = left_side_exists.at(it.first);
    auto rse = right_side_exists.at(it.first);

    if (!lse || !rse) {
      continue;
    }

    auto lsf = left_side_face.at(it.first);
    auto rsf = right_side_face.at(it.first);

    auto segment = it.first.second;

    auto &left_ste = face_to_local_segment_to_eid.at(lsf);
    auto &right_ste = face_to_local_segment_to_eid.at(rsf);

    auto left_eid = left_ste.at(segment);
    auto right_eid = right_ste.at(segment);

    auto left_ef = face_to_local_EF.at(lsf).row(left_eid);
    auto right_ef = face_to_local_EF.at(rsf).row(right_eid);

    auto left_face = left_ef(0) == -1 ? left_ef(1) : left_ef(0);
    auto right_face = right_ef(0) == -1 ? right_ef(1) : right_ef(0);

    auto left_face_id = face_to_tessellation_unified.at(lsf)(left_face);
    auto right_face_id = face_to_tessellation_unified.at(rsf)(right_face);

    glueing_pairs.push_back({left_face_id, right_face_id});
  }

  return glueing_pairs;
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXi>
glue_tessellations(const Eigen::MatrixXd &V_unified,
                   const Eigen::MatrixXi &F_unified,
                   const std::vector<std::pair<int, int>> &glueing_pairs) {
  Eigen::MatrixXi F_unified_new = F_unified;
  for (auto pair : glueing_pairs) {
    int f0 = pair.first;
    int f1 = pair.second;

    // for each vertex in f0, find the closest vertex in f1
    for (int i = 0; i < 3; i++) {
      double min_dist = std::numeric_limits<double>::infinity();
      int min_v1 = -1;
      for (int j = 0; j < 3; j++) {
        int v0 = F_unified(f0, i);
        int v1 = F_unified(f1, j);

        auto v0_pos = V_unified.row(v0);
        auto v1_pos = V_unified.row(v1);

        auto dist = (v0_pos - v1_pos).norm();
        if (dist < min_dist) {
          min_dist = dist;
          min_v1 = v1;
        }
      }

      if (min_v1 != -1) {
        throw std::runtime_error("Glueing pairs are not valid");
      }

      F_unified_new(f0, i) = min_v1;
    }
  }

  // remove redundant vertices
  Eigen::MatrixXd V_unified_new;
  Eigen::VectorXi I, J;
  igl::remove_unreferenced(V_unified, F_unified_new, V_unified_new,
                           F_unified_new, I, J);

  return {V_unified_new, F_unified_new, J};
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, std::vector<Point_2>>
glue_tessellations(
    const std::map<std::pair<int, Face_const_handle>, tessellation_t>
        &face_to_tessellation,
    const std::map<std::pair<int, Segment_2>, bool, IntSegmentPairComparator>
        &left_side_exists,
    const std::map<std::pair<int, Segment_2>, bool, IntSegmentPairComparator>
        &right_side_exists,
    const std::map<std::pair<int, Segment_2>, std::pair<int, Face_const_handle>,
                   IntSegmentPairComparator> &left_side_face,
    const std::map<std::pair<int, Segment_2>, std::pair<int, Face_const_handle>,
                   IntSegmentPairComparator> &right_side_face
) {
  // 1. unify everything into a single mesh
  auto [face_to_tessellation_unified, V_unified, F_unified, V_unified_to_point] = unify_tessellations(face_to_tessellation);

  // 2. get the list of faces to glue on the global tessellation
  auto glueing_pairs = get_glueing_pairs(
      face_to_tessellation, face_to_tessellation_unified, left_side_exists,
      right_side_exists, left_side_face, right_side_face);

  // 3. glue the tessellations
  auto [V_glued, F_glued, V_glued_to_unified] =
      glue_tessellations(V_unified, F_unified, glueing_pairs);

  std::vector<Point_2> V_glued_to_point(V_glued.rows());
  for (int i = 0; i < V_glued.rows(); i++) {
    V_glued_to_point[i] = V_unified_to_point[V_glued_to_unified(i)];
  }

  // 3. return the result
  return {V_glued, F_glued, V_glued_to_point};
}
} // namespace utils