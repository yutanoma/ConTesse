#include "assign_qi.h"

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

bool is_identical_segment(const Segment_2 &s1, const Segment_2 &s2) {
    return (s1.source() == s2.source() && s1.target() == s2.target()) ||
           (s1.source() == s2.target() && s1.target() == s2.source());
}

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


// tells if the given point is a crossing.
// if it is, it also returns the qi for each side of the crossing
std::tuple<bool, std::vector<std::pair<Segment_2, int>>> is_crossing(
    const Arrangement_with_history_2 &arr,
    const std::map<Segment_2, int, Segment2Comparator> &qi,
    const std::map<Point_2, Segment_2> &upper_casing_edges,
    const std::map<Point_2, Segment_2> &lower_casing_edges,
    const std::map<Point_2, Point_2> &next,
    Vertex_const_iterator &vertex
) {
    std::vector<std::pair<Segment_2, int>> returning_qi;
    if (vertex->degree() != 4) {
        return {false, returning_qi};
    }

    auto vertex_point = vertex->point();

    std::vector<Segment_2> upper_new_edges, lower_new_edges;
    Segment_2 upper_old_edge, lower_old_edge;

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
            upper_old_edge = segment;
        } else if (is_identical_segment(segment, lower_casing_edges.at(vertex->point()))) {
            lower_new_edges.push_back(it->curve());
            lower_old_edge = segment;
        }

        ++it;
    } while (it != vertex->incident_halfedges());

    if (upper_new_edges.size() != 2 || lower_new_edges.size() != 2) {
        std::cerr << "Error: upper_new_edges.size(): " << upper_new_edges.size()
                  << " and lower_new_edges.size(): " << lower_new_edges.size()
                  << " edges" << std::endl;
        throw std::runtime_error("Error: upper_new_edges.size(): " + std::to_string(upper_new_edges.size()) + " and lower_new_edges.size(): " + std::to_string(lower_new_edges.size()) + " edges");
    }

    Point_2 upper_1_pt_new = upper_new_edges[0].source() == vertex_point ? upper_new_edges[0].target() : upper_new_edges[0].source();
    Point_2 upper_2_pt_new = upper_new_edges[1].source() == vertex_point ? upper_new_edges[1].target() : upper_new_edges[1].source();
    Point_2 lower_1_pt_new = lower_new_edges[0].source() == vertex_point ? lower_new_edges[0].target() : lower_new_edges[0].source();
    Point_2 lower_2_pt_new = lower_new_edges[1].source() == vertex_point ? lower_new_edges[1].target() : lower_new_edges[1].source();

    Segment_2 upper_1_new = upper_new_edges[0];
    Segment_2 upper_2_new = upper_new_edges[1];
    Segment_2 lower_1_new = lower_new_edges[0];
    Segment_2 lower_2_new = lower_new_edges[1];

    // these are the prev and next points of the old edges
    Point_2 upper_old_prev, upper_old_next;
    Point_2 lower_old_prev, lower_old_next;

    if (next.at(upper_old_edge.source()) == upper_old_edge.target()) {
        upper_old_prev = upper_old_edge.source();
        upper_old_next = upper_old_edge.target();
    } else if (next.at(upper_old_edge.target()) == upper_old_edge.source()) {
        upper_old_prev = upper_old_edge.target();
        upper_old_next = upper_old_edge.source();
    } else {
        std::cerr << "Error: upper_old_edge.source(): " << upper_old_edge.source() << " and upper_old_edge.target(): " << upper_old_edge.target() << std::endl;
        throw std::runtime_error("Error: upper_old_edge.source() and upper_old_edge.target() are not the same");
    }
    
    if (next.at(lower_old_edge.source()) == lower_old_edge.target()) {
        lower_old_prev = lower_old_edge.source();
        lower_old_next = lower_old_edge.target();
    } else if (next.at(lower_old_edge.target()) == lower_old_edge.source()) {
        lower_old_prev = lower_old_edge.target();
        lower_old_next = lower_old_edge.source();
    } else {
        std::cerr << "Error: lower_old_edge.source(): " << lower_old_edge.source() << " and lower_old_edge.target(): " << lower_old_edge.target() << std::endl;
        throw std::runtime_error("Error: lower_old_edge.source() and lower_old_edge.target() are not the same");
    }

    // the closer one to _old_prev is the _new_prev
    Point_2 upper_new_prev, upper_new_next;
    
    if (
        CGAL::squared_distance(upper_1_pt_new, upper_old_prev) <
        CGAL::squared_distance(upper_2_pt_new, upper_old_prev)
    ) {
        upper_new_prev = upper_1_pt_new;
        upper_new_next = upper_2_pt_new;
    } else {
        upper_new_prev = upper_2_pt_new;
        upper_new_next = upper_1_pt_new;
    }

    // Given the vector from upper_prev to upper_next,
    // the one on the right side of the vector is the right side of upper_1 and upper_2
    Eigen::Vector2d upper_vec(
        CGAL::to_double(upper_new_next.x() - upper_new_prev.x()), 
        CGAL::to_double(upper_new_next.y() - upper_new_prev.y())
    );

    Eigen::Vector2d lower_1_vec(
        CGAL::to_double(lower_1_pt_new.x() - upper_new_prev.x()), 
        CGAL::to_double(lower_1_pt_new.y() - upper_new_prev.y())
    );
    Eigen::Vector2d lower_2_vec(
        CGAL::to_double(lower_2_pt_new.x() - upper_new_prev.x()), 
        CGAL::to_double(lower_2_pt_new.y() - upper_new_prev.y())
    );

    // Compute 2D cross product (z-component of 3D cross product)
    double cross_1 = upper_vec.x() * lower_1_vec.y() - upper_vec.y() * lower_1_vec.x();
    double cross_2 = upper_vec.x() * lower_2_vec.y() - upper_vec.y() * lower_2_vec.x();

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

    if (
        qi.at(lower_visible) != -1 &&
        qi.at(lower_visible) >= upper_qi && // upper casing edges must have higher qi than lower casing edges
        upper_qi != -1 // upper casing edges must have a qi
    ) {
        int qi_base = qi.at(lower_visible);
        returning_qi.push_back({lower_visible, qi_base});
        returning_qi.push_back({lower_occluded, qi_base + 1});
    } else if (
        qi.at(lower_occluded) != -1 &&
        qi.at(lower_occluded) > 0 && // qi_base - 1 > 0 must hold
        qi.at(lower_occluded) >= upper_qi && // upper casing edges must have higher qi than lower casing edges
        upper_qi != -1 // upper casing edges must have a qi
    ) {
        int qi_base = qi.at(lower_occluded);
        returning_qi.push_back({lower_visible, qi_base - 1});
        returning_qi.push_back({lower_occluded, qi_base});
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
    const std::map<Point_2, Segment_2> &upper_casing_edges,
    const std::map<Point_2, Segment_2> &lower_casing_edges,
    const std::map<Point_2, Point_2> &next,
    std::map<Segment_2, int, Segment2Comparator> &qi
) {
    // the queue of segments to be propagated
    SegmentQueue q;

    // 1. find out all vertices where adjacent edges are yet to be assigned qi
    std::map<Point_2, Vertex_const_iterator> point_to_vertex;
    for (auto it = arr.vertices_begin(); it != arr.vertices_end(); it++) {
        auto vertex = it;

        point_to_vertex[vertex->point()] = vertex;

        auto [crossing_found, qi_crossing] = is_crossing(arr, qi, upper_casing_edges, lower_casing_edges, next, vertex);

        if (crossing_found) {
            for (auto [segment, qi_value] : qi_crossing) {
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
                auto [crossing_found, qi_crossing] = is_crossing(arr, qi, upper_casing_edges, lower_casing_edges, next, vertex);

                if (crossing_found) {
                    for (auto [segment, qi_value] : qi_crossing) {
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

                // {
                //     Point_2 dbg(-0.0813, 0.3326);
                //     if (CGAL::squared_distance(dbg, current_vertex->point()) < 1e-6) {
                //         std::cout << "current_vertex->point(): " << current_vertex->point() << std::endl;
                //         std::cout << "current_segment: " << current_segment << std::endl;
                //         std::cout << "current_qi: " << current_qi << std::endl;

                //         std::cout << "qi_crossing: " << qi_crossing.size() << std::endl;
                //         for (auto [segment, qi_value] : qi_crossing) {
                //             std::cout << "  segment: " << segment << ", qi_value: " << qi_value << std::endl;
                //             std::cout << "      qi.at(segment): " << qi.at(segment) << std::endl;
                //         }
                //         std::cout << std::endl;
                //     }
                // }

                break;
            }
        }
    }
}

std::tuple<bool, std::vector<Point_2>> check_qi_mismatch(
    const Arrangement_with_history_2 &arr,
    const std::map<Point_2, bool> &point_is_singularity,
    const std::map<Segment_2, bool, Segment2Comparator> &segment_is_convex,
    // given an intersection point, which edges are the upper and lower casing edges?
    const std::map<Point_2, Segment_2> &upper_casing_edges,
    const std::map<Point_2, Segment_2> &lower_casing_edges,
    // given a point, which point is the next point in the circular order?
    // *exclude* the self-intersection point in this list
    const std::map<Point_2, Point_2> &next,
    const std::map<Segment_2, int, Segment2Comparator> &qi
) {
    bool is_valid = true;
    std::vector<Point_2> qi_mismatch_positions;
    
    // check all vertices in the arrangement
    for (auto it = arr.vertices_begin(); it != arr.vertices_end(); it++) {
        auto vertex = it;

        auto [crossing_found, qi_crossing] = is_crossing(arr, qi, upper_casing_edges, lower_casing_edges, next, vertex);
        if (crossing_found) {
            if (qi_crossing.size() != 4) {
                // the qi must be assigned to all 4 sides respecting the casing information
                // qi_crossing.size() != 4 means that some sides are not assigned a qi that conforms to the casing information
                is_valid = false;
                qi_mismatch_positions.push_back(vertex->point());
            } else {
                for (auto [segment, qi_value] : qi_crossing) {
                    if (qi.at(segment) != qi_value) {
                        is_valid = false;
                        qi_mismatch_positions.push_back(vertex->point());
                        break;
                    }
                }
            }
        }
        
        auto [singularity_found, qi_singularity] = is_singularity(arr, point_is_singularity, segment_is_convex, qi, vertex);
        if (singularity_found) {
            for (auto [segment, qi_value] : qi_singularity) {
                if (qi.at(segment) != qi_value) {
                    is_valid = false;
                    qi_mismatch_positions.push_back(vertex->point());
                    break;
                }
            }
        }
    }

    return {is_valid, qi_mismatch_positions};
}

// @return: whether the arrangement is valid, qi for each segment, and positions where qi mismatch occurs
std::tuple<bool, std::map<Segment_2, int, Segment2Comparator>, std::vector<Point_2>> validity_check(
    Arrangement_with_history_2 &arr,
    std::map<Point_2, bool> &point_is_singularity,
    // convex/concave labeling per segment. This is based on the "old" segments on the planar map
    std::map<Segment_2, bool, Segment2Comparator> &segment_is_convex,
    // given an intersection point, which "old" edges are the upper and lower casing edges?
    std::map<Point_2, Segment_2> &upper_casing_edges,
    std::map<Point_2, Segment_2> &lower_casing_edges,
    // given a point, which point is the next point in the circular order?
    // *exclude* the self-intersection point in this list
    std::map<Point_2, Point_2> &next
) {
    // store the propagated qi
    // this is based on the "new" segments on the planar map
    std::map<Segment_2, int, Segment2Comparator> qi;
    for (auto it = arr.edges_begin(); it != arr.edges_end(); it++) {
        Segment_2 segment = it->curve();
        qi[segment] = -1;
    }

    // 1. propagate the qi from the silhouette
    assign_silhouette_qi(arr, qi);
    propagate_qi(arr, point_is_singularity, segment_is_convex, upper_casing_edges, lower_casing_edges, next, qi);

    // // next, for each independent hole, try to assign the qi
    // // the qi must be equipped with the depth of the hole
    // std::map<Segment_2, int, Segment2Comparator> depth;
    // for (auto it = arr.edges_begin(); it != arr.edges_end(); it++) {
    //     Segment_2 segment = it->curve();
    //     depth[segment] = 0;
    //     // the qi is depth + qi
    // }

    // assign_hole_qi(arr, qi, depth);

    // // propagate the qi from the holes
    // propagate_qi(arr, qi);

    // // determine the d for each hole
    // determine_hole_d(arr, qi, depth);

    // finally, check the validity for each vertex / singularity
    auto [is_valid, qi_mismatch_positions] = check_qi_mismatch(arr, point_is_singularity, segment_is_convex, upper_casing_edges, lower_casing_edges, next, qi);

    // // also check if the winding number constraints are met
    // bool is_valid_winding_number = check_winding_number_constraints(arr, qi, depth);
    // is_valid = is_valid && is_valid_winding_number;

    return {is_valid, qi, qi_mismatch_positions};
}
}
