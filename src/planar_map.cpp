#include "planar_map.h"
#include "segment_contours.h"
#include <map>
#include <iostream>
#include <cmath>
#include <algorithm>

#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_with_history_2.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace utils {
// Helper function to check if a point lies on a segment
bool point_on_segment(const Point_2& p, const Point_2& seg_start, const Point_2& seg_end) {
    // Check if point p lies on segment from seg_start to seg_end
    // First check collinearity, then check if p is between seg_start and seg_end
    if (CGAL::collinear(seg_start, p, seg_end)) {
        // Check if p is between seg_start and seg_end using coordinate comparison
        auto min_x = std::min(CGAL::to_double(seg_start.x()), CGAL::to_double(seg_end.x()));
        auto max_x = std::max(CGAL::to_double(seg_start.x()), CGAL::to_double(seg_end.x()));
        auto min_y = std::min(CGAL::to_double(seg_start.y()), CGAL::to_double(seg_end.y()));
        auto max_y = std::max(CGAL::to_double(seg_start.y()), CGAL::to_double(seg_end.y()));
        
        auto px = CGAL::to_double(p.x());
        auto py = CGAL::to_double(p.y());
        
        return (px >= min_x && px <= max_x && py >= min_y && py <= max_y);
    }
    return false;
}

// Helper function to find which original segment a halfedge belongs to
std::tuple<int, int, int, bool> find_original_segment_info(
    const Point_2& source, const Point_2& target,
    const std::vector<std::vector<Point_2>>& original_segments) {
    
    // Check each original segment to see if this halfedge lies along it
    for (int loop_id = 0; loop_id < original_segments.size(); ++loop_id) {
        const auto& loop_points = original_segments[loop_id];
        for (int edge_id = 0; edge_id < loop_points.size(); ++edge_id) {
            int next_id = (edge_id + 1) % loop_points.size();
            const Point_2& seg_start = loop_points[edge_id];
            const Point_2& seg_end = loop_points[next_id];
            
            // Check if both source and target lie on this original segment
            if (point_on_segment(source, seg_start, seg_end) && 
                point_on_segment(target, seg_start, seg_end)) {
                
                // Determine direction: is source->target same direction as seg_start->seg_end?
                bool is_forward = true;
                
                // Use parameter values to determine direction along the segment
                auto seg = Segment_2(seg_start, seg_end);
                auto source_param = CGAL::to_double((source.x() - seg_start.x()) * (seg_end.x() - seg_start.x()) + 
                                                  (source.y() - seg_start.y()) * (seg_end.y() - seg_start.y())) /
                                  CGAL::to_double((seg_end.x() - seg_start.x()) * (seg_end.x() - seg_start.x()) + 
                                                (seg_end.y() - seg_start.y()) * (seg_end.y() - seg_start.y()));
                auto target_param = CGAL::to_double((target.x() - seg_start.x()) * (seg_end.x() - seg_start.x()) + 
                                                  (target.y() - seg_start.y()) * (seg_end.y() - seg_start.y())) /
                                  CGAL::to_double((seg_end.x() - seg_start.x()) * (seg_end.x() - seg_start.x()) + 
                                                (seg_end.y() - seg_start.y()) * (seg_end.y() - seg_start.y()));
                
                is_forward = (target_param > source_param);
                
                return std::make_tuple(loop_id * 1000 + edge_id, loop_id, edge_id, is_forward);
            }
        }
    }
    
    return std::make_tuple(-1, -1, -1, true); // Not found
}

// Helper function to compute winding number of a point with respect to a polygon
// todo: avoid traversing all the points in the polygon
double compute_winding_number(const Point_2& query_point, const std::vector<Point_2>& polygon) {
    double winding_number = 0.0;
    int n = polygon.size();
    
    for (int i = 0; i < n; ++i) {
        int next_i = (i + 1) % n;
        
        Point_2 p1 = polygon[i];
        Point_2 p2 = polygon[next_i];
        
        // Convert to double coordinates for computation
        double x1 = CGAL::to_double(p1.x());
        double y1 = CGAL::to_double(p1.y());
        double x2 = CGAL::to_double(p2.x());
        double y2 = CGAL::to_double(p2.y());
        double qx = CGAL::to_double(query_point.x());
        double qy = CGAL::to_double(query_point.y());
        
        // Translate so query point is at origin
        x1 -= qx; y1 -= qy;
        x2 -= qx; y2 -= qy;
        
        // Compute the signed angle contribution
        // Using the atan2 method for robustness
        double angle1 = std::atan2(y1, x1);
        double angle2 = std::atan2(y2, x2);
        
        // Compute angle difference, handling wraparound
        double diff = angle2 - angle1;
        
        // Normalize to [-π, π]
        while (diff > M_PI) diff -= 2 * M_PI;
        while (diff < -M_PI) diff += 2 * M_PI;
        
        winding_number += diff;
    }
    
    // Convert from radians to winding number (divide by 2π)
    return winding_number / (2.0 * M_PI);
}

// Helper function to test if a point is inside a polygon using winding number
bool is_point_inside_polygon(const Point_2& point, const std::vector<Point_2>& polygon) {
    if (polygon.size() < 3) return false;
    
    // Compute winding number of the polygon boundary with respect to the test point
    // For a counterclockwise-oriented polygon boundary:
    // - winding number ≈ +1 means point is inside
    // - winding number ≈ 0 means point is outside
    double winding = compute_winding_number(point, polygon);
    
    // Use a tolerance for numerical robustness
    return std::abs(winding - 1.0) < 0.5; // Should be close to +1 if inside
}

// Helper function to find a representative point inside a bounded face
Point_2 find_face_representative_point(const Arrangement_2::Face_const_handle& face) {
    if (face->is_unbounded()) {
        // For unbounded face, return a point far outside
        return Point_2(-1000, -1000);
    }
    
    // Collect boundary points in order
    std::vector<Point_2> boundary_points;
    auto first = face->outer_ccb();
    auto curr = first;
    do {
        boundary_points.push_back(curr->source()->point());
        ++curr;
    } while (curr != first);
    
    if (boundary_points.size() < 3) {
        return Point_2(0, 0); // Degenerate case
    }
    
    // Find an interior point using the "ear" method:
    // For any three consecutive vertices that form a valid triangle,
    // the centroid of that triangle should be inside the face
    for (int i = 0; i < boundary_points.size(); ++i) {
        int i1 = i;
        int i2 = (i + 1) % boundary_points.size();
        int i3 = (i + 2) % boundary_points.size();
        
        Point_2 p1 = boundary_points[i1];
        Point_2 p2 = boundary_points[i2];
        Point_2 p3 = boundary_points[i3];
        
        // Check if this forms a "left turn" (interior angle)
        // Using CGAL's orientation test
        if (CGAL::orientation(p1, p2, p3) == CGAL::LEFT_TURN) {
            // Compute triangle centroid
            double cx = (CGAL::to_double(p1.x()) + CGAL::to_double(p2.x()) + CGAL::to_double(p3.x())) / 3.0;
            double cy = (CGAL::to_double(p1.y()) + CGAL::to_double(p2.y()) + CGAL::to_double(p3.y())) / 3.0;
            
            Point_2 candidate(cx, cy);
            
            // Verify this point is inside the polygon using a simple ray casting test
            if (is_point_inside_polygon(candidate, boundary_points)) {
                return candidate;
            }
        }
    }
    
    // Fallback: use centroid but move it slightly inward
    double sum_x = 0, sum_y = 0;
    for (const auto& p : boundary_points) {
        sum_x += CGAL::to_double(p.x());
        sum_y += CGAL::to_double(p.y());
    }
    
    Point_2 centroid(sum_x / boundary_points.size(), sum_y / boundary_points.size());
    
    // If centroid is inside, use it; otherwise try to move it inward
    if (is_point_inside_polygon(centroid, boundary_points)) {
        return centroid;
    }
    
    // Last resort: find any interior point by averaging with boundary midpoints
    for (int i = 0; i < boundary_points.size(); ++i) {
        int next_i = (i + 1) % boundary_points.size();
        Point_2 p1 = boundary_points[i];
        Point_2 p2 = boundary_points[next_i];
        
        // Midpoint of edge
        double mx = (CGAL::to_double(p1.x()) + CGAL::to_double(p2.x())) / 2.0;
        double my = (CGAL::to_double(p1.y()) + CGAL::to_double(p2.y())) / 2.0;
        
        // Average with centroid and move slightly inward
        double tx = (mx + sum_x / boundary_points.size()) / 2.0;
        double ty = (my + sum_y / boundary_points.size()) / 2.0;
        
        Point_2 candidate(tx, ty);
        if (is_point_inside_polygon(candidate, boundary_points)) {
            return candidate;
        }
    }
    
    return centroid; // Fallback
}

// Helper function to calculate t parameter for a point along an original segment
double t_parameter(const Point_2& point, const Point_2& seg_start, const Point_2& seg_end) {
    double dx_seg = CGAL::to_double(seg_end.x() - seg_start.x());
    double dy_seg = CGAL::to_double(seg_end.y() - seg_start.y());
    double length = std::sqrt(dx_seg * dx_seg + dy_seg * dy_seg);
    
    if (length < 1e-10) return 0.0;
    
    double dx_point = CGAL::to_double(point.x() - seg_start.x());
    double dy_point = CGAL::to_double(point.y() - seg_start.y());
    double distance = std::sqrt(dx_point * dx_point + dy_point * dy_point);
    
    double t = distance / length;
    return std::max(0.0, std::min(1.0, t));
}

// New function that works with arrangement history to get correct t values
std::vector<Intersection> planar_map_intersections(
    const Arrangement_with_history_2 &arr_with_history
) {
    std::vector<Intersection> intersections;

    // Iterate through arrangement vertices to find intersection points
    for (auto vit = arr_with_history.vertices_begin(); vit != arr_with_history.vertices_end(); ++vit) {
        if (vit->is_isolated() || vit->degree() <= 2) {
            continue;  // Not an intersection vertex
        }

        Point_2 vertex_point = vit->point();

        auto first_he = vit->incident_halfedges();
        auto curr_he = first_he;

        std::set<Segment_2, Segment2Comparator> original_segments;

        do {
            // For each incident halfedge, find which original curves it belongs to
            if (arr_with_history.number_of_originating_curves(curr_he) != 1) {
                std::cerr << "Error: curr_he has " << arr_with_history.number_of_originating_curves(curr_he)
                          << " originating curves" << std::endl;
                continue;
            }

            auto ocit = arr_with_history.originating_curves_begin(curr_he);
            Segment_2 curve_handle = *ocit;

            original_segments.insert(curve_handle);

            ++curr_he;
        } while (curr_he != first_he);
        
        if (original_segments.size() != 2) {
            std::cerr << "Error: curr_he has " << original_segments.size()
                      << " originating segments" << std::endl;
            continue;
        }

        Intersection intersection;
        intersection.segment1 = *original_segments.begin();
        intersection.segment2 = *std::next(original_segments.begin());

        Point_2 s1_start = intersection.segment1.source();
        Point_2 s1_end = intersection.segment1.target();
        intersection.t1 = t_parameter(vertex_point, s1_start, s1_end);

        intersection.p1_1 = s1_start;
        intersection.p1_2 = s1_end;

        Point_2 s2_start = intersection.segment2.source();
        Point_2 s2_end = intersection.segment2.target();
        intersection.t2 = t_parameter(vertex_point, s2_start, s2_end);

        intersection.p2_1 = s2_start;
        intersection.p2_2 = s2_end;

        intersection.intersection_point = vertex_point;

        intersections.push_back(intersection);
    }
    
    return intersections;
}

std::tuple<
    Arrangement_with_history_2,
    std::map<Point_2, std::pair<int, int>>,
    std::map<Segment_2, std::pair<int, int>, Segment2Comparator>
> planar_map(
    const std::vector<
        std::pair<Eigen::MatrixXd, Eigen::MatrixXi>
    > &contours_2d
) {
    // Create CGAL Arrangement_2 object
    Arrangement_with_history_2 arr;
    
    // Store original segments for later reference
    std::vector<std::vector<Point_2>> original_segments;
    
    // First pass: collect all unique points and store original segments
    std::vector<Segment_2> segments_to_insert;

    // these maps are used to keep track of the original segments
    // in the pair, the first int is the curve index
    // the second int is the node index
    std::map<Point_2, std::pair<int, int>> point_to_node;

    // these maps are used to keep track of the original segments
    // in the pair, the first int is the curve index
    // the second int is the segment index
    std::map<Segment_2, std::pair<int, int>, Segment2Comparator> segment_to_segment;

    for (int i = 0; i < contours_2d.size(); ++i) {
        const auto& contour = contours_2d[i];
        const auto& N = std::get<0>(contour);
        const auto& S = std::get<1>(contour);

        std::unordered_map<int, Point_2> point_map;

        for (int j = 0; j < N.rows(); ++j) {
            auto p = Point_2(N(j, 0), N(j, 1));
            point_map.insert({j, p});
            point_to_node.insert({p, {i, j}});
        }

        for (int j = 0; j < S.rows(); ++j) {
            int start_idx = S(j, 0);
            int end_idx = S(j, 1);

            Point_2 p1 = point_map[start_idx];
            Point_2 p2 = point_map[end_idx];

            Segment_2 segment(p1, p2);
            segments_to_insert.push_back(segment);
            segment_to_segment.insert({segment, {i, j}});
        }
    }

    // Insert all segments into the arrangement
    CGAL::insert(arr, segments_to_insert.begin(), segments_to_insert.end());

    return std::make_tuple(arr, point_to_node, segment_to_segment);
}
}
