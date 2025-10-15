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

using SegmentQueue = std::priority_queue<QueueEntry, std::vector<QueueEntry>, CompareByQi>;

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
std::tuple<bool, std::vector<std::pair<Segment_2, int>>> is_crossing(
    const Arrangement_with_history_2 &arr,
    const std::map<Segment_2, int, Segment2Comparator> &qi,
    const std::map<Point_2, Segment_2> &upper_casing_edges,
    const std::map<Point_2, Segment_2> &lower_casing_edges,
    const std::map<Segment_2, int, Segment2Comparator> &segment_orientation,
    const std::map<Segment_2, bool, Segment2Comparator> &is_cut,
    Vertex_const_iterator &vertex
) {
    std::vector<std::pair<Segment_2, int>> returning_qi;
    if (vertex->degree() != 4) {
        return {false, returning_qi};
    }

    auto vertex_point = vertex->point();

    std::vector<Segment_2> upper_new_edges, lower_new_edges;

    auto it = vertex->incident_halfedges();
    do {
        auto segment_new = it->curve();
        
        // this must have an old edge
        if (arr.number_of_originating_curves(it) != 1) {
            std::cerr << "Error: it has " << arr.number_of_originating_curves(it)
                      << " originating curves" << std::endl;
            throw std::runtime_error("Error: it has " + std::to_string(arr.number_of_originating_curves(it)) + " originating curves");
        }

        auto ocit = arr.originating_curves_begin(it);
        Segment_2 segment = *ocit;

        if (is_identical_segment(segment, upper_casing_edges.at(vertex->point()))) {
            upper_new_edges.push_back(it->curve());
        } else if (is_identical_segment(segment, lower_casing_edges.at(vertex->point()))) {
            lower_new_edges.push_back(it->curve());
        }

        ++it;
    } while (it != vertex->incident_halfedges());

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

    bool is_upper_1_positive = segment_orientation.at(upper_1_new) == 1;
    bool is_upper_2_positive = segment_orientation.at(upper_2_new) == 1;
    bool is_lower_1_positive = segment_orientation.at(lower_1_new) == 1;
    bool is_lower_2_positive = segment_orientation.at(lower_2_new) == 1;

    bool upper_1_prev_crossing = is_upper_1_positive ? upper_1_new.source() == vertex_point : upper_1_new.target() == vertex_point;
    bool upper_2_prev_crossing = is_upper_2_positive ? upper_2_new.source() == vertex_point : upper_2_new.target() == vertex_point;
    bool lower_1_prev_crossing = is_lower_1_positive ? lower_1_new.source() == vertex_point : lower_1_new.target() == vertex_point;
    bool lower_2_prev_crossing = is_lower_2_positive
                                     ? lower_2_new.source() == vertex_point
                                     : lower_2_new.target() == vertex_point;

    if ((upper_1_prev_crossing && upper_2_prev_crossing) || (!upper_1_prev_crossing && !upper_2_prev_crossing)) {
      std::cerr << "Error: upper_1_prev_crossing and upper_2_prev_crossing are "
                   "both true or both false"
                << std::endl;
      std::cerr << "upper_1_prev_crossing: " << upper_1_prev_crossing << ", upper_2_prev_crossing: " << upper_2_prev_crossing << std::endl;
      throw std::runtime_error("Error: upper_1_prev_crossing and upper_2_prev_crossing are both true");
    }
    if ((lower_1_prev_crossing && lower_2_prev_crossing) || (!lower_1_prev_crossing && !lower_2_prev_crossing)) {
      std::cerr << "Error: lower_1_prev_crossing and lower_2_prev_crossing are "
                   "both true or both false"
                << std::endl;
      std::cerr << "lower_1_prev_crossing: " << lower_1_prev_crossing << ", lower_2_prev_crossing: " << lower_2_prev_crossing << std::endl;
      throw std::runtime_error("Error: lower_1_prev_crossing and lower_2_prev_crossing are both true");
    }

    Segment_2 upper_prev_edge = upper_1_prev_crossing ? upper_2_new : upper_1_new;
    Segment_2 upper_next_edge = upper_1_prev_crossing ? upper_1_new : upper_2_new;
    Segment_2 lower_prev_edge = lower_1_prev_crossing ? lower_2_new : lower_1_new;
    Segment_2 lower_next_edge = lower_1_prev_crossing ? lower_1_new : lower_2_new;

    // these are the points that are adjacent to the vertex in question
    Point_2 upper_prev_pt = upper_prev_edge.source() == vertex_point ? upper_prev_edge.target() : upper_prev_edge.source();
    Point_2 upper_next_pt = upper_next_edge.source() == vertex_point ? upper_next_edge.target() : upper_next_edge.source();
    Point_2 lower_prev_pt = lower_prev_edge.source() == vertex_point ? lower_prev_edge.target() : lower_prev_edge.source();
    Point_2 lower_next_pt = lower_next_edge.source() == vertex_point ? lower_next_edge.target() : lower_next_edge.source();

    // Given the vector from upper_prev to upper_next,
    // the one on the right side of the vector is the right side of upper_1 and upper_2
    Eigen::Vector2d upper_vec(
        CGAL::to_double(upper_next_pt.x() - upper_prev_pt.x()), 
        CGAL::to_double(upper_next_pt.y() - upper_prev_pt.y())
    );

    // a vector from vertex to lower_prev
    Eigen::Vector2d lower_prev_vec(
        CGAL::to_double(lower_prev_pt.x() - vertex_point.x()), 
        CGAL::to_double(lower_prev_pt.y() - vertex_point.y())
    );
    Eigen::Vector2d lower_next_vec(
        CGAL::to_double(lower_next_pt.x() - vertex_point.x()), 
        CGAL::to_double(lower_next_pt.y() - vertex_point.y())
    );

    // Compute 2D cross product (z-component of 3D cross product)
    double cross_1 = upper_vec.x() * lower_prev_vec.y() - upper_vec.y() * lower_prev_vec.x();
    double cross_2 = upper_vec.x() * lower_next_vec.y() - upper_vec.y() * lower_next_vec.x();

    Segment_2 lower_occluded, lower_visible;

    if (cross_1 < 0 && cross_2 > 0) {
        // 1 is on the right side and 2 is on the left side
        // right side means that the qi is lower / it is visible
        lower_occluded = lower_2_new;
        lower_visible = lower_1_new;
    } else if (cross_1 > 0 && cross_2 < 0) {
        // 2 is on the right side and 1 is on the left side
        // right side means that the qi is lower / it is visible
        lower_occluded = lower_1_new;
        lower_visible = lower_2_new;
    } else {
        std::cout << "cross_1: " << cross_1 << ", cross_2: " << cross_2 << std::endl;
        throw std::runtime_error("Error: both are on the same side");
    }

    // among upper_1, upper_2, lower_occluded, lower_visible,
    // find if any have qi
    int upper_qi = -1;
    if (qi.at(upper_1_new) != -1) {
        upper_qi = qi.at(upper_1_new);
        returning_qi.push_back({upper_1_new, upper_qi});
        returning_qi.push_back({upper_2_new, upper_qi});
    } else if (qi.at(upper_2_new) != -1) {
        upper_qi = qi.at(upper_2_new);
        returning_qi.push_back({upper_1_new, upper_qi});
        returning_qi.push_back({upper_2_new, upper_qi});
    }

    int lower_qi = -1;

    if (is_cut.at(upper_1_new) && is_cut.at(upper_2_new)) {
      // if the upper is a cut, nothing changes on all the points
      lower_qi = qi.at(lower_visible);

      if (lower_qi != -1) {
        returning_qi.push_back({lower_visible, lower_qi});
        returning_qi.push_back({lower_occluded, lower_qi});
      }
    } else if (is_cut.at(lower_1_new) && is_cut.at(lower_2_new)) {
      // if the lower is a cut, only the cut gets smaller on the left side of the upper
      lower_qi = qi.at(lower_visible);

      if (lower_qi != -1) {
        returning_qi.push_back({lower_visible, lower_qi});
        returning_qi.push_back({lower_occluded, lower_qi + 1});
      }
    } else {
      lower_qi = qi.at(lower_visible);

      if (lower_qi != -1 && lower_qi >= upper_qi) {
        returning_qi.push_back({lower_visible, lower_qi});
        returning_qi.push_back({lower_occluded, lower_qi + 1});
      } 
    }

    return {true, returning_qi};
}

std::tuple<bool, std::vector<std::pair<Segment_2, int>>>
is_threeway(const Arrangement_with_history_2 &arr,
            const std::map<Segment_2, int, Segment2Comparator> &qi,
            const std::map<Segment_2, bool, Segment2Comparator> &is_cut,
            Vertex_const_iterator &vertex) {
  // a threeway has one cut and two non-cut edges
  std::vector<Segment_2> cut_segments;
  std::vector<Segment_2> non_cut_segments;
  for (auto it = vertex->incident_halfedges(); it != vertex->incident_halfedges(); it++) {
    if (is_cut.at(it->curve())) {
      cut_segments.push_back(it->curve());
    } else {
      non_cut_segments.push_back(it->curve());
    }
  }
  
  if (cut_segments.size() != 1 || non_cut_segments.size() != 2) {
    return {false, {}};
  }

  // on a threeway, the qi must be assigned as q, q, q
  int qi_value = -1;
  if (qi.at(cut_segments[0]) != -1) {
    qi_value = qi.at(cut_segments[0]);
  } else if (qi.at(non_cut_segments[0]) != -1) {
    qi_value = qi.at(non_cut_segments[0]);
  } else if (qi.at(non_cut_segments[1]) != -1) {
    qi_value = qi.at(non_cut_segments[1]);
  }

  std::vector<std::pair<Segment_2, int>> returning_qi;

  if (qi_value != -1) {
    returning_qi.push_back({cut_segments[0], qi_value});
    returning_qi.push_back({non_cut_segments[0], qi_value});
    returning_qi.push_back({non_cut_segments[1], qi_value});
  }

  return {true, returning_qi};
}

// tells if the given point is a singularity.
// if it is, it also returns the qi for each side of the singularity
std::tuple<bool, std::vector<std::pair<Segment_2, int>>> is_singularity(
    const Arrangement_with_history_2 &arr,
    const std::map<Point_2, bool> &point_is_singularity,
    const std::map<Segment_2, bool, Segment2Comparator> &segment_is_convex,
    const std::map<Segment_2, int, Segment2Comparator> &qi,
    Vertex_const_iterator &vertex
) {
    std::vector<std::pair<Segment_2, int>> returning_qi;

    if (vertex->degree() != 2) {
        return {false, returning_qi};
    }

    if (!point_is_singularity.at(vertex->point())) {
        return {false, returning_qi};
    }

    Segment_2 occluded, visible;

    if (vertex->degree() != 2) {
        std::cerr << "Error: vertex has " << vertex->degree() << " edges" << std::endl;
        throw std::runtime_error("Error: vertex has " + std::to_string(vertex->degree()) + " edges");
    }

    auto it = vertex->incident_halfedges();
    do {
        // we must use the originating segments for querying segment_is_convex
        
        if (arr.number_of_originating_curves(it) != 1) {
            std::cerr << "Error: it has " << arr.number_of_originating_curves(it)
                      << " originating curves" << std::endl;
            throw std::runtime_error("Error: it has " + std::to_string(arr.number_of_originating_curves(it)) + " originating curves");
        }

        auto ocit = arr.originating_curves_begin(it);

        Segment_2 segment = *ocit;

        if (segment_is_convex.at(segment)) {
            visible = it->curve();
        } else {
            occluded = it->curve();
        }

        ++it;
    } while (it != vertex->incident_halfedges());

    if (qi.at(visible) != -1) {
        returning_qi.push_back({visible, qi.at(visible)});
        returning_qi.push_back({occluded, qi.at(visible) + 1});
    } else if (qi.at(occluded) != -1 && qi.at(occluded) > 0) {
        returning_qi.push_back({occluded, qi.at(occluded)});
        returning_qi.push_back({visible, qi.at(occluded) - 1});
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
    const std::map<Point_2, Segment_2> &upper_casing_edges,
    const std::map<Point_2, Segment_2> &lower_casing_edges,
    const std::map<Segment_2, int, Segment2Comparator> &segment_orientation,
    std::map<Segment_2, int, Segment2Comparator> &qi
) {
    // the queue of segments to be propagated
    SegmentQueue q;

    // 1. find out all vertices where adjacent edges are yet to be assigned qi
    std::map<Point_2, Vertex_const_iterator> point_to_vertex;
    for (auto it = arr.vertices_begin(); it != arr.vertices_end(); it++) {
        auto vertex = it;

        point_to_vertex[vertex->point()] = vertex;

        auto [crossing_found, qi_crossing] = is_crossing(arr, qi, upper_casing_edges, lower_casing_edges, segment_orientation, segment_is_cut, vertex);

        if (crossing_found) {
            for (auto [segment, qi_value] : qi_crossing) {
                if (qi.at(segment) == -1) {
                    q.push({qi_value, segment});
                }
            }
        }

        auto [threeway_found, qi_threeway] = is_threeway(arr, qi, segment_is_cut, vertex);

        if (threeway_found) {
            for (auto [segment, qi_value] : qi_threeway) {
                if (qi.at(segment) == -1) {
                    q.push({qi_value, segment});
                }
            }
        }

        auto [singularity_found, qi_singularity] = is_singularity(arr, point_is_singularity, segment_is_convex, qi, vertex);

        if (singularity_found) {
            for (auto [segment, qi_value] : qi_singularity) {
                if (qi.at(segment) == -1) {
                    q.push({qi_value, segment});
                }
            }

        }
    }
    
    while (!q.empty()) {
        Segment_2 current_segment = q.top().segment;
        int current_qi = q.top().qi;
        q.pop();

        if (qi.at(current_segment) != -1) {
            continue;
        }

        qi[current_segment] = current_qi;

        auto v1 = point_to_vertex[current_segment.source()];
        auto v2 = point_to_vertex[current_segment.target()];

        Vertex_const_iterator current_vertex;
        bool is_both_termini = false;

        if (is_terminus(arr, point_is_singularity, v1) && is_terminus(arr, point_is_singularity, v2)) {
            current_vertex = v1;
            is_both_termini = true;
        } else if (!is_terminus(arr, point_is_singularity, v1)) {
            current_vertex = v1;
        } else if (!is_terminus(arr, point_is_singularity, v2)) {
            current_vertex = v2;
        } else {
            throw std::runtime_error("Error: neither are termini");
        }

        // 2. for each assigned side, propagate the qi until it hits another ambiguous point
        // if that ambiguous point has no -1 qi sides, then stop propagating
        while (!is_terminus(arr, point_is_singularity, current_vertex)) {
            auto prev_current_segment = current_segment;

            auto it = current_vertex->incident_halfedges();
            do {
                Segment_2 segment = it->curve();

                if (!is_identical_segment(current_segment, segment)) {
                    current_segment = segment;
                    break;
                }

                ++it;
            } while (it != current_vertex->incident_halfedges());

            auto next_v = current_vertex->point() == current_segment.source() ? current_segment.target() : current_segment.source();
            current_vertex = point_to_vertex[next_v];

            qi[current_segment] = current_qi;
        }

        std::vector<Vertex_const_iterator> vertices_to_check;
        if (is_both_termini) {
            vertices_to_check.push_back(v1);
            vertices_to_check.push_back(v2);
        } else {
            vertices_to_check.push_back(current_vertex);
        }

        for (auto vertex : vertices_to_check) {
            if (is_terminus(arr, point_is_singularity, vertex)) {
                auto [crossing_found, qi_crossing] = is_crossing(arr, qi, upper_casing_edges, lower_casing_edges, segment_orientation, segment_is_cut, vertex);

                if (crossing_found) {
                    for (auto [segment, qi_value] : qi_crossing) {
                        if (qi.at(segment) == -1) {
                            q.push({qi_value, segment});
                        }
                    }
                }

                auto [threeway_found, qi_threeway] = is_threeway(arr, qi, segment_is_cut, vertex);

                if (threeway_found) {
                    for (auto [segment, qi_value] : qi_threeway) {
                        if (qi.at(segment) == -1) {
                            q.push({qi_value, segment});
                        }
                    }
                }

                auto [singularity_found, qi_singularity] = is_singularity(arr, point_is_singularity, segment_is_convex, qi, vertex);
        
                if (singularity_found) {
                    for (auto [segment, qi_value] : qi_singularity) {
                        if (qi.at(segment) == -1) {
                            q.push({qi_value, segment});
                        }
                    }
                }

                break;
            }
        }
    }
}

// @return: whether the arrangement is valid, qi for each segment, and positions where qi mismatch occurs
std::tuple<bool, std::map<Segment_2, int, Segment2Comparator>, std::vector<Point_2>> validity_check(
    Arrangement_with_history_2 &arr,
    std::map<Point_2, bool> &point_is_singularity,
    // convex/concave labeling per segment. This is based on the "old" segments on the planar map
    std::map<Segment_2, bool, Segment2Comparator> &segment_is_convex,
    // the segment is a cut on the triangulation
    std::map<Segment_2, bool, Segment2Comparator> &segment_is_cut,
    // given an intersection point, which "old" edges are the upper and lower casing edges?
    std::map<Point_2, Segment_2> &upper_casing_edges,
    std::map<Point_2, Segment_2> &lower_casing_edges,
    // given a segment, does segment.source() -> segment.target() matches the canonical orientation?
    std::map<Segment_2, int, Segment2Comparator> &segment_orientation
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
                 upper_casing_edges, lower_casing_edges, segment_orientation, qi);

    // 2. propagate the qi from the edges where the winding number on the right side face is the lowest
    assign_lowest_wn_qi(arr, wn, segment_orientation, qi);
    propagate_qi(arr, point_is_singularity, segment_is_convex, segment_is_cut,
                 upper_casing_edges, lower_casing_edges, segment_orientation, qi);

    // finally, check the validity for each vertex / singularity
    std::map<Segment_2, bool, Segment2Comparator> is_segment_invalid;
    for (auto it = arr.edges_begin(); it != arr.edges_end(); it++) {
        Segment_2 segment = it->curve();
        is_segment_invalid[segment] = false;
    }

    auto [is_valid, qi_mismatch_positions] = check_qi_mismatch(arr, segment_orientation, qi, is_segment_invalid);

    return {is_valid, qi, qi_mismatch_positions};
}
} // namespace utils
