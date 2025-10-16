#include "validity_check.h"
#include "planar_map.h"

#include <Eigen/Geometry>

namespace utils {
std::tuple <
    std::map<Segment_2, Face_const_handle, Segment2Comparator>, // left face
    std::map<Segment_2, Face_const_handle, Segment2Comparator>  // right face
> segment_to_face(const Arrangement_with_history_2 &arr,
  const std::map<Segment_2, int, Segment2Comparator> &segment_orientation) {
  
  std::map<Segment_2, Face_const_handle, Segment2Comparator> segment_to_left_face;
  std::map<Segment_2, Face_const_handle, Segment2Comparator> segment_to_right_face;
  
  // Iterate through all edges in the arrangement
  for (auto it = arr.edges_begin(); it != arr.edges_end(); it++) {
    auto segment = it->curve();
    auto orientation_it = segment_orientation.find(segment);
    
    if (orientation_it == segment_orientation.end()) {
      // If segment not found, try the reverse direction
      auto reverse_segment = Segment_2(segment.target(), segment.source());
      orientation_it = segment_orientation.find(reverse_segment);
    }
    
    if (orientation_it != segment_orientation.end()) {
      int orientation = orientation_it->second;
      
      // Get the two faces adjacent to this edge
      auto left_face = it->face();
      auto right_face = it->twin()->face();
      
      // Determine if halfedge direction matches segment direction
      bool directions_match = (it->source()->point() == segment.source() && it->target()->point() == segment.target());
      
      Face_const_handle actual_left_face, actual_right_face;
      
      if (orientation == 1) {
        // Left side gets +1, right side gets -1
        if (directions_match) {
          // Current face is on the left, twin face is on the right
          actual_left_face = left_face;
          actual_right_face = right_face;
        } else {
          // Current face is on the right, twin face is on the left
          actual_left_face = right_face;
          actual_right_face = left_face;
        }
      } else if (orientation == -1) {
        // Right side gets +1, left side gets 0 (flipped from above)
        if (directions_match) {
          // Current face is on the right, twin face is on the left
          actual_left_face = right_face;
          actual_right_face = left_face;
        } else {
          // Current face is on the left, twin face is on the right
          actual_left_face = left_face;
          actual_right_face = right_face;
        }
      }
      
      // Store the mapping
      segment_to_left_face[segment] = actual_left_face;
      segment_to_right_face[segment] = actual_right_face;
    }
  }
  
  return {segment_to_left_face, segment_to_right_face};
}

std::map<Face_const_handle, int> face_wn(const Arrangement_with_history_2 &arr,
                                         const std::map<Segment_2, int, Segment2Comparator> &segment_orientation
) {
  std::map<Face_const_handle, int> wn;
  
  // Step 1: Initialize with unbounded face (winding number = 0)
  auto unbounded_face = arr.unbounded_face();
  wn[unbounded_face] = 0;
  
  // Step 2: Use BFS to propagate winding numbers
  std::queue<Face_const_handle> face_queue;
  face_queue.push(unbounded_face);

  std::unordered_map<Face_const_handle, bool> visited;
  for (auto it = arr.faces_begin(); it != arr.faces_end(); it++) {
    visited[it] = false;
  }

  auto [segment_to_left_face, segment_to_right_face] = segment_to_face(arr, segment_orientation);
  
  size_t processed = 0;
  while (face_queue.size() != 0) {
    auto current_face = face_queue.front();
    face_queue.pop();

    if (visited[current_face]) {
      continue;
    }

    visited[current_face] = true;

    // Process all incident edges (both outer boundary and holes)
    auto process_halfedge = [&](const auto &he) {
      auto twin_face = he->twin()->face();

      if (visited[twin_face]) {
        return;
      }

      auto segment = he->curve();
      auto left_face = segment_to_left_face.at(segment);
      auto right_face = segment_to_right_face.at(segment);

      if (left_face == current_face) {
        wn[right_face] = wn[current_face] - 1;
      } else if (right_face == current_face) {
        wn[left_face] = wn[current_face] + 1;
      } else {
        throw std::runtime_error("Error: neither are on the same side");
      }

      face_queue.push(twin_face);
    };

    // Process outer boundary
    if (current_face->has_outer_ccb()) {
      auto outer_ccb = current_face->outer_ccb();
      if (outer_ccb != CGAL::Arrangement_2<Traits_2>::Ccb_halfedge_circulator()) {
        auto he = outer_ccb;
        do {
          process_halfedge(he);
          ++he;
        } while (he != outer_ccb);
      }
    }

    // Process holes
    for (auto hole_it = current_face->holes_begin(); hole_it != current_face->holes_end(); ++hole_it) {
      auto hole_ccb = *hole_it;
      auto he = hole_ccb;
      do {
        process_halfedge(he);
        ++he;
      } while (he != hole_ccb);
    }
    
    processed++;
  }
  
  return wn;
}

// check if the qi is consistent with the arrangement
// it checks if
// (1) the qi matches the topological constraints of an image-space intersection, and
// (2) the qi matches the winding number constraints
// no need for the convexity or casing information
std::tuple<bool, std::vector<Point_2>>
check_qi_mismatch(const Arrangement_with_history_2 &arr,
                  // if the actual orientation is segement.source() -> segement.target(),
                  // then the orientation is 1. otherwise, it is -1.
                  const std::map<Segment_2, int, Segment2Comparator> &segment_orientation,
                  const std::map<Segment_2, int, Segment2Comparator> &qi,
                  // invalid segments will be ignored for this check
                  const std::map<Segment_2, bool, Segment2Comparator> &is_segment_invalid
) {
  bool is_valid = true;
  std::vector<Point_2> qi_mismatch_positions;

  // 1. check all vertices in the arrangement
  for (auto it = arr.vertices_begin(); it != arr.vertices_end(); it++) {
    auto vertex = it;

    // 1-1. at an image space intersection, the qi must be assigned as q, q, p, p+1, where q <= p
    if (vertex->degree() == 4) {
      bool is_image_space_intersection = true;
      for (auto it = vertex->incident_halfedges(); it != vertex->incident_halfedges(); it++) {
        auto segment = it->curve();

        if (is_segment_invalid.at(segment)) {
          is_image_space_intersection = false;
          break;
        }
      }

      if (!is_image_space_intersection) {
        continue;
      }

      std::vector<int> qi_values;
      for (auto it = vertex->incident_halfedges(); it != vertex->incident_halfedges(); it++) {
        auto segment = it->curve();
        qi_values.push_back(qi.at(segment));
      }

      int qi_max = *std::max_element(qi_values.begin(), qi_values.end());
      int qi_min = *std::min_element(qi_values.begin(), qi_values.end());

      // q has to be qi_min and p+1 has to be qi_max
      int q = qi_min;
      int p = qi_max - 1;
      int p_plus_1 = qi_max;

      int q_count = 0;
      int p_count = 0;
      int p_plus_1_count = 0;
      for (auto qi_value : qi_values) {
        if (qi_value == q) {
          q_count++;
        } else if (qi_value == p) {
          p_count++;
        } else if (qi_value == p_plus_1) {
          p_plus_1_count++;
        }
      }

      if (q == p) {
        // there has to be three q's and one p+1
        q_count = q_count + p_count;
        p_count = 0;
        if (q_count != 3 || p_plus_1_count != 1) {
          is_valid = false;
          qi_mismatch_positions.push_back(vertex->point());
        }
      } else {
        // there has to be two q's and one p+1 and one p
        if (q_count != 2 || p_count != 1 || p_plus_1_count != 1) {
          is_valid = false;
          qi_mismatch_positions.push_back(vertex->point());
        }
      }
    }

    // 1-2. at a singularity, the qi must be assigned as q, q+1
    if (vertex->degree() == 2) {
      bool is_singularity = true;
      std::vector<int> qi_values;
      for (auto it = vertex->incident_halfedges(); it != vertex->incident_halfedges(); it++) {
        auto segment = it->curve();
        if (is_segment_invalid.at(segment)) {
          is_singularity = false;
          break;
        }

        qi_values.push_back(qi.at(segment));
      }

      if (!is_singularity || qi_values.size() != 2) {
        continue;
      }

      if (qi_values[0] != qi_values[1]) {
        int diff = qi_values[1] - qi_values[0];
        if (std::abs(diff) != 1) {
          is_valid = false;
          qi_mismatch_positions.push_back(vertex->point());
        }
      }
    }
  }

  // 2. check all edges in the arrangement
  // if the qi is between zero and the winding number of the right side face
  auto [segment_to_left_face, segment_to_right_face] = segment_to_face(arr, segment_orientation);
  auto wn = face_wn(arr, segment_orientation);
  for (auto it = arr.edges_begin(); it != arr.edges_end(); it++) {
    auto segment = it->curve();

    auto right_face = segment_to_right_face.at(segment);
    auto right_wn = wn.at(right_face);

    auto qi_value = qi.at(segment);

    // the qi_value must be between zero and the winding number of the right side face
    if (qi_value < 0 || qi_value > right_wn) {
      is_valid = false;
      qi_mismatch_positions.push_back(it->source()->point());
      qi_mismatch_positions.push_back(it->target()->point());
    }
  }

  return {is_valid, qi_mismatch_positions};
}
} // namespace utils
