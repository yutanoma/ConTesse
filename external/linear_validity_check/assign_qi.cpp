#include "assign_qi.h"

#include "validity_check.h"
#include "planar_map.h"

#include <Eigen/Geometry>

namespace utils {
typedef CGAL::Arrangement_with_history_2<Traits_2>::Vertex_handle Vertex_handle;
typedef CGAL::Arrangement_with_history_2<Traits_2>::Vertex_const_iterator Vertex_const_iterator;

struct QueueEntry {
    int qi;
    Segment_2 segment;
};

struct CompareByQi {
bool operator()(const QueueEntry& a, const QueueEntry& b) const noexcept {
    // return true if 'a' has lower priority than 'b'
    return a.qi > b.qi; // min-heap on qi: smaller qi comes first
}
};

using SegmentQueue = std::queue<QueueEntry>;

void assign_silhouette_qi(
    Arrangement_with_history_2 &arr,
    std::map<Segment_2, int, Segment2Comparator> &qi
) {
    auto unbounded_face = arr.unbounded_face();
    // The unbounded face has no outer CCB. Its boundary components are stored as holes.
    for (auto hit = unbounded_face->holes_begin(); hit != unbounded_face->holes_end(); ++hit) {
        auto ccb = *hit; // CCB circulator for a hole (boundary of a bounded component)
        auto curr = ccb;

        if (curr == nullptr) {
            std::cerr << "Error: curr is nullptr" << std::endl;
            throw std::runtime_error("Error: curr is nullptr");
        }

        do {
            Segment_2 segment = curr->curve();
            qi[segment] = 0;

            ++curr;  // move to next half-edge in the CCB
        } while (curr != ccb);
    }
}

void assign_lowest_wn_qi(
    const Arrangement_with_history_2 &arr,
    const std::map<Face_const_handle, int> &wn,
    const std::map<Segment_2, int, Segment2Comparator> &segment_orientation,
    std::map<Segment_2, int, Segment2Comparator> &qi
) {
  std::map<Segment_2, bool, Segment2Comparator> visited;
  std::queue<Segment_2> q;

  for (auto it = arr.edges_begin(); it != arr.edges_end(); it++) {
    if (qi.at(it->curve()) != -1) {
      visited[it->curve()] = true;
    } else {
      visited[it->curve()] = false;
      q.push(it->curve());
    }
  }

  std::map<Point_2, std::vector<Segment_2>> point_to_segments;
  for (auto it = arr.edges_begin(); it != arr.edges_end(); it++) {
    point_to_segments[it->source()->point()].push_back(it->curve());
    point_to_segments[it->target()->point()].push_back(it->curve());
  }

  auto [left_face, right_face] = segment_to_face(arr, segment_orientation);

  while (!q.empty()) {
    auto segment = q.front();
    q.pop();

    if (visited.at(segment)) {
      continue;
    }

    int smallest_wn = wn.at(right_face.at(segment));
    Segment_2 smallest_wn_segment = segment;

    // iterate through all the connected segments
    std::queue<Segment_2> q_connected;
    q_connected.push(segment);

    while (!q_connected.empty()) {
      auto segment_connected = q_connected.front();
      q_connected.pop();

      visited[segment_connected] = true;

      int right_wn = wn.at(right_face.at(segment_connected));
      if (right_wn < smallest_wn) {
        smallest_wn = right_wn;
        smallest_wn_segment = segment_connected;
      }

      // add the incident segments
      std::vector<Point_2> incident_points = {segment_connected.source(),
                                              segment_connected.target()};
      for (auto incident_point : incident_points) {
        auto segments = point_to_segments.at(incident_point);

        for (auto segment_connected : segments) {
          if (!visited.at(segment_connected)) {
            q_connected.push(segment_connected);
          }
        }
      }
    }
  }
}


// tells if the given point is a crossing.
// if it is, it also returns the qi for each side of the crossing
std::tuple<bool, std::vector<std::pair<Halfedge_const_handle, int>>> is_crossing(
    const Arrangement_with_history_2 &arr,
    const std::map<Segment_2, int, Segment2Comparator> &qi,
    const std::map<Point_2, std::vector<Segment_2>> &upper_casing_edges,
    const std::map<Point_2, std::vector<Segment_2>> &lower_casing_edges,
    const std::map<Segment_2, int, Segment2Comparator> &segment_orientation,
    const std::map<Segment_2, bool, Segment2Comparator> &is_cut,
    Vertex_const_iterator &vertex
) {
    std::vector<std::pair<Halfedge_const_handle, int>> returning_qi;
    if (vertex->degree() != 4) {
        return {false, returning_qi};
    }

    auto vertex_point = vertex->point();

    if (upper_casing_edges.find(vertex_point) == upper_casing_edges.end() ||
        lower_casing_edges.find(vertex_point) == lower_casing_edges.end()) {
        std::cerr << "Error: upper_casing_edges or lower_casing_edges not found for vertex: " << vertex_point << std::endl;
        throw std::runtime_error("Error: upper_casing_edges or lower_casing_edges not found for vertex");
    }

    std::vector<Segment_2> upper_new_edges = upper_casing_edges.at(vertex_point);
    std::vector<Segment_2> lower_new_edges = lower_casing_edges.at(vertex_point);

    if (upper_new_edges.size() != 2 || lower_new_edges.size() != 2) {
        std::cerr << "Error: upper_new_edges.size(): " << upper_new_edges.size()
                  << " and lower_new_edges.size(): " << lower_new_edges.size()
                  << " edges" << std::endl;
        throw std::runtime_error("Error: upper_new_edges.size(): " + std::to_string(upper_new_edges.size()) + " and lower_new_edges.size(): " + std::to_string(lower_new_edges.size()) + " edges");
    }

    Segment_2 upper_1_new = upper_new_edges[0];
    Segment_2 upper_2_new = upper_new_edges[1];
    Segment_2 lower_1_new = lower_new_edges[0];
    Segment_2 lower_2_new = lower_new_edges[1];

    Halfedge_const_handle upper_1_edge, upper_2_edge, lower_1_edge, lower_2_edge;
    for (auto edge_it = vertex->incident_halfedges(); edge_it != vertex->incident_halfedges(); edge_it++) {
      auto segment = edge_it->curve();
      if (is_identical_segment(segment, upper_1_new)) {
        upper_1_edge = edge_it;
      } else if (is_identical_segment(segment, upper_2_new)) {
        upper_2_edge = edge_it;
      } else if (is_identical_segment(segment, lower_1_new)) {
        lower_1_edge = edge_it;
      } else if (is_identical_segment(segment, lower_2_new)) {
        lower_2_edge = edge_it;
      }
    }

    bool is_upper_1_positive = segment_orientation.at(upper_1_new) == 1;
    
    // Given the vector upper_1, check if lower_1 and lower_2 are on the right and left sides of the vector
    Eigen::Vector2d upper_vec;

    if (is_upper_1_positive) {
      upper_vec = Eigen::Vector2d(
        CGAL::to_double(upper_1_new.target().x() - upper_1_new.source().x()),
        CGAL::to_double(upper_1_new.target().y() - upper_1_new.source().y())
      );
    } else {
      upper_vec = Eigen::Vector2d(
        CGAL::to_double(upper_1_new.source().x() - upper_1_new.target().x()),
        CGAL::to_double(upper_1_new.source().y() - upper_1_new.target().y())
      );
    }

    // a vector from vertex to lower_1/lower_2
    Point_2 lower_1_pt = is_identical_point(lower_1_new.source(), vertex_point) ? lower_1_new.target() : lower_1_new.source();
    Point_2 lower_2_pt = is_identical_point(lower_2_new.source(), vertex_point) ? lower_2_new.target() : lower_2_new.source();

    Eigen::Vector2d lower_1_vec(
        CGAL::to_double(lower_1_pt.x() - vertex_point.x()), 
        CGAL::to_double(lower_1_pt.y() - vertex_point.y())
    );
    Eigen::Vector2d lower_2_vec(
        CGAL::to_double(lower_2_pt.x() - vertex_point.x()), 
        CGAL::to_double(lower_2_pt.y() - vertex_point.y())
    );

    // Compute 2D cross product (z-component of 3D cross product)
    // if this is positive, then lower is on the left side of the vector
    double cross_1 = upper_vec.x() * lower_1_vec.y() - upper_vec.y() * lower_1_vec.x();
    double cross_2 = upper_vec.x() * lower_2_vec.y() - upper_vec.y() * lower_2_vec.x();

    Halfedge_const_handle lower_occluded_edge, lower_visible_edge;

    if (cross_1 < 0 && cross_2 > 0) {
        // 1 is on the right side and 2 is on the left side
        // right side means that the qi is lower / it is visible
        lower_occluded_edge = lower_2_edge;
        lower_visible_edge = lower_1_edge;
    } else if (cross_1 > 0 && cross_2 < 0) {
        // 2 is on the right side and 1 is on the left side
        // right side means that the qi is lower / it is visible
        lower_occluded_edge = lower_1_edge;
        lower_visible_edge = lower_2_edge;
    } else {
      std::cout << "  cross_1: " << cross_1 << ", cross_2: " << cross_2
                << std::endl;
      std::cout << "  upper_vec: " << upper_vec.transpose() << std::endl;
      std::cout << "  lower_1_vec: " << lower_1_vec.transpose() << std::endl;
      std::cout << "  lower_2_vec: " << lower_2_vec.transpose() << std::endl;
        throw std::runtime_error("Error: both are on the same side");
    }

    // among upper_1, upper_2, lower_occluded, lower_visible,
    // find if any have qi
    int upper_qi = -1;
    if (qi.at(upper_1_new) != -1) {
        upper_qi = qi.at(upper_1_new);
        returning_qi.push_back({upper_1_edge, upper_qi});
        returning_qi.push_back({upper_2_edge, upper_qi});
    } else if (qi.at(upper_2_new) != -1) {
        upper_qi = qi.at(upper_2_new);
        returning_qi.push_back({upper_1_edge, upper_qi});
        returning_qi.push_back({upper_2_edge, upper_qi});
    }

    int lower_qi = -1;

    if (is_cut.at(upper_1_new) && is_cut.at(upper_2_new)) {
      // if the upper is a cut, nothing changes on all the lower points
      if (qi.at(lower_visible_edge->curve()) != -1) {
        lower_qi = qi.at(lower_visible_edge->curve());
      } else if (qi.at(lower_occluded_edge->curve()) != -1) {
        lower_qi = qi.at(lower_occluded_edge->curve());
      }

      if (lower_qi != -1) {
        returning_qi.push_back({lower_visible_edge, lower_qi});
        returning_qi.push_back({lower_occluded_edge, lower_qi});
      }
    } else {
      // otherwise, adopt the usual rules
      if (qi.at(lower_visible_edge->curve()) != -1) {
        lower_qi = qi.at(lower_visible_edge->curve()) - 1;
        returning_qi.push_back({lower_visible_edge, lower_qi});
        returning_qi.push_back({lower_occluded_edge, lower_qi + 1});
      } else if (qi.at(lower_occluded_edge->curve()) != -1) {
        lower_qi = qi.at(lower_occluded_edge->curve());
        returning_qi.push_back({lower_visible_edge, lower_qi - 1});
        returning_qi.push_back({lower_occluded_edge, lower_qi});
      }
    }

    return {true, returning_qi};
}

std::tuple<bool, std::vector<std::pair<Halfedge_const_handle, int>>>
is_threeway(const Arrangement_with_history_2 &arr,
            const std::map<Segment_2, int, Segment2Comparator> &qi,
            const std::map<Segment_2, bool, Segment2Comparator> &is_cut,
            Vertex_const_iterator &vertex) {
  if (vertex->degree() != 3) {
    return {false, {}};
  }
  
  // a threeway has one cut and two non-cut edges
  int cut_count = 0, non_cut_count = 0;
  std::vector<Halfedge_const_handle> all_edges;
  auto it = vertex->incident_halfedges();
  do {
    if (is_cut.at(it->curve())) {
      cut_count++;
    } else {
      non_cut_count++;
    }
    all_edges.push_back(it);
    ++it;
  } while (it != vertex->incident_halfedges());

  if (cut_count != 1 && cut_count != 3) {
    return {false, {}};
  }

  // on a threeway, the qi must be assigned as q, q, q
  std::vector<std::pair<Halfedge_const_handle, int>> returning_qi;
  int qi_value = -1;
  for (auto edge_it : all_edges) {
    auto segment = edge_it->curve();
    if (qi.at(segment) != -1) {
      qi_value = qi.at(segment);
      break;
    }
  }

  if (qi_value != -1) {
    returning_qi.push_back({all_edges[0], qi_value});
    returning_qi.push_back({all_edges[1], qi_value});
    returning_qi.push_back({all_edges[2], qi_value});
  }

  return {true, returning_qi};
}

// tells if the given point is a singularity.
// if it is, it also returns the qi for each side of the singularity
std::tuple<bool, std::vector<std::pair<Halfedge_const_handle, int>>> is_singularity(
    const Arrangement_with_history_2 &arr,
    const std::map<Point_2, bool> &point_is_singularity,
    const std::map<Segment_2, bool, Segment2Comparator> &segment_is_convex,
    const std::map<Segment_2, int, Segment2Comparator> &qi,
    const bool is_back_facing,
    Vertex_const_iterator &vertex
) {
    std::vector<std::pair<Halfedge_const_handle, int>> returning_qi;

    if (vertex->degree() != 2) {
        return {false, returning_qi};
    }

    if (!point_is_singularity.at(vertex->point())) {
        return {false, returning_qi};
    }

    Halfedge_const_handle occluded, visible;

    if (vertex->degree() != 2) {
        std::cerr << "Error: vertex has " << vertex->degree() << " edges" << std::endl;
        throw std::runtime_error("Error: vertex has " + std::to_string(vertex->degree()) + " edges");
    }

    auto it = vertex->incident_halfedges();
    do {
        // we must use the originating segments for querying segment_is_convex
        auto segment = it->curve();

        if (segment_is_convex.at(segment)) {
            visible = it;
        } else {
            occluded = it;
        }

        ++it;
    } while (it != vertex->incident_halfedges());

    // on a back facing patch, the visible side is concave and the occluded side is convex
    if (is_back_facing) {
      std::swap(visible, occluded);
      // std::cout << "  visible: " << visible << ", occluded: " << occluded << std::endl;
    }

    if (qi.at(visible->curve()) != -1) {
        returning_qi.push_back({visible, qi.at(visible->curve())});
        returning_qi.push_back({occluded, qi.at(visible->curve()) + 1});
    } else if (qi.at(occluded->curve()) != -1 && qi.at(occluded->curve()) > 0) {
        returning_qi.push_back({occluded, qi.at(occluded->curve())});
        returning_qi.push_back({visible, qi.at(occluded->curve()) - 1});
    } else {
        // nothing can be said
    }

    return {true, returning_qi};
}

bool is_terminus(
    const Arrangement_with_history_2 &arr,
    const std::map<Point_2, bool> &point_is_singularity,
    Vertex_const_iterator &vertex
) {
    return vertex->degree() != 2 || point_is_singularity.at(vertex->point());
}

void propagate_qi(
    const Arrangement_with_history_2 &arr,
    const std::map<Point_2, bool> &point_is_singularity,
    const std::map<Segment_2, bool, Segment2Comparator> &segment_is_convex,
    const std::map<Segment_2, bool, Segment2Comparator> &segment_is_cut,
    const std::map<Point_2, std::vector<Segment_2>> &upper_casing_edges,
    const std::map<Point_2, std::vector<Segment_2>> &lower_casing_edges,
    const std::map<Segment_2, int, Segment2Comparator> &segment_orientation,
    const bool is_back_facing,
    std::map<Segment_2, int, Segment2Comparator> &qi
) {
    // the queue of segments to be propagated
    std::queue<std::pair<int, Halfedge_const_handle>> q;

    // add the segments with qi != -1
    for (auto it = arr.edges_begin(); it != arr.edges_end(); it++) {
      auto segment = it->curve();
      if (qi.at(segment) != -1) {
        q.push({qi.at(segment), it});
      }
    }

    while (!q.empty()) {
      auto [current_qi, edge_it] = q.front();
      q.pop();

      auto segment = edge_it->curve();

      // add its neighboring segments
      std::vector<Arrangement_with_history_2::Vertex_const_handle>
          adjacent_vertices = {edge_it->source(), edge_it->target()};
      for (auto vertex : adjacent_vertices) {
        bool _is_terminus = is_terminus(arr, point_is_singularity, vertex);
        if (_is_terminus) {
          auto [_is_sing, _sing_qi] = is_singularity(arr, point_is_singularity, segment_is_convex, qi, is_back_facing, vertex);
          if (_is_sing) {
            for (auto [edge_it, qi_value] : _sing_qi) {
              if (qi_value != -1 && qi.at(edge_it->curve()) == -1) {
                  qi[edge_it->curve()] = qi_value;
                  q.push({qi_value, edge_it});
              }
            }
          }

          auto [_is_threeway, _three_qi] =
              is_threeway(arr, qi, segment_is_cut, vertex);
          if (_is_threeway) {
            for (auto [edge_it, qi_value] : _three_qi) {
              std::cout << "  threeway: " << edge_it->curve() << ", qi_value: " << qi_value << std::endl;
              if (qi_value != -1 && qi.at(edge_it->curve()) == -1) {
                  qi[edge_it->curve()] = qi_value;
                  q.push({qi_value, edge_it});
              }
            }
          }

          auto [_is_crossing, _crossing_qi] =
              is_crossing(arr, qi, upper_casing_edges, lower_casing_edges, segment_orientation, segment_is_cut, vertex);
          if (_is_crossing) {
            for (auto [edge_it, qi_value] : _crossing_qi) {
              if (qi_value != -1 && qi.at(edge_it->curve()) == -1) {
                  qi[edge_it->curve()] = qi_value;
                  q.push({qi_value, edge_it});
              }
            }
          }
        } else {
          auto e1 = vertex->incident_halfedges();
          auto e2 = vertex->incident_halfedges()++;

          if (is_identical_segment(e1->curve(), segment) && qi.at(e2->curve()) == -1) {
            qi[e2->curve()] = current_qi;
            q.push({current_qi, e2});
          } else if (is_identical_segment(e2->curve(), segment) && qi.at(e1->curve()) == -1) {
            qi[e1->curve()] = current_qi;
            q.push({current_qi, e1});
          }
        }
      }
    }
}

// @return: whether the arrangement is valid, qi for each segment, and positions where qi mismatch occurs
std::tuple<bool, std::map<Segment_2, int, Segment2Comparator>, std::vector<Point_2>> validity_check(
    Arrangement_with_history_2 &arr,
    std::map<Point_2, bool> &point_is_singularity,
    // convex/concave labeling per segment. This is based on the "new" segments on the planar map
    std::map<Segment_2, bool, Segment2Comparator> &segment_is_convex,
    // the segment is a cut on the triangulation
    std::map<Segment_2, bool, Segment2Comparator> &segment_is_cut,
    // given an intersection point, which "new" edges are the upper and lower casing edges?
    std::map<Point_2, std::vector<Segment_2>> &upper_casing_edges,
    std::map<Point_2, std::vector<Segment_2>> &lower_casing_edges,
    // given a "new" segment, does segment.source() -> segment.target() matches the canonical orientation?
    std::map<Segment_2, int, Segment2Comparator> &segment_orientation,
    bool is_back_facing
) {
    // store the propagated qi
    // this is based on the "new" segments on the planar map
    std::map<Segment_2, int, Segment2Comparator> qi;
    for (auto it = arr.edges_begin(); it != arr.edges_end(); it++) {
        Segment_2 segment = it->curve();
        qi[segment] = -1;
    }

    // 0. assign the winding number to each face
    std::map<Face_const_handle, int> wn = face_wn(arr, segment_orientation);

    // 1. propagate the qi from the silhouette
    assign_silhouette_qi(arr, qi);
    propagate_qi(arr, point_is_singularity, segment_is_convex, segment_is_cut,
                 upper_casing_edges, lower_casing_edges, segment_orientation, is_back_facing, qi);

    // 2. propagate the qi from the edges where the winding number on the right side face is the lowest
    assign_lowest_wn_qi(arr, wn, segment_orientation, qi);
    propagate_qi(arr, point_is_singularity, segment_is_convex, segment_is_cut,
                 upper_casing_edges, lower_casing_edges, segment_orientation, is_back_facing, qi);

    // finally, check the validity for each vertex / singularity
    std::map<Segment_2, bool, Segment2Comparator> is_segment_invalid;
    for (auto it = arr.edges_begin(); it != arr.edges_end(); it++) {
        Segment_2 segment = it->curve();
        if (segment_is_cut.at(segment)) {
          is_segment_invalid[segment] = true;
        } else {
          is_segment_invalid[segment] = false;
        }
    }

    auto [is_valid, qi_mismatch_positions] = check_qi_mismatch(arr, segment_orientation, qi, is_segment_invalid);

    return {is_valid, qi, qi_mismatch_positions};
}
} // namespace utils
