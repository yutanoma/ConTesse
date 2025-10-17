#include "planar_map_to_svg.h"
#include "planar_map.h"
#include "validity_check.h"

#include <sstream>
#include <iomanip>
#include <unordered_set>
#include <limits>
#include <cmath>
#include <algorithm>
#include <set>

namespace utils {
static std::string color_for_qi(int q) {
    static const std::vector<std::string> palette = {
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
        "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
        "#bcbd22", "#17becf"
    };


    if (q < 0) {
        return "#444444";
    }

    int idx = (q >= 0 ? q : -q) % static_cast<int>(palette.size());
    return palette[idx];
}

static bool same_segment(const Segment_2 &a, const Segment_2 &b) {
    return (a.source() == b.source() && a.target() == b.target()) ||
           (a.source() == b.target() && a.target() == b.source());
}

template <typename T>
static bool map_lookup_segment_bi_oriented(
    const std::map<Segment_2, T, Segment2Comparator> &m,
    const Segment_2 &s,
    T &out_val
) {
    auto it = m.find(s);
    if (it != m.end()) { out_val = it->second; return true; }
    Segment_2 rs(s.target(), s.source());
    it = m.find(rs);
    if (it != m.end()) { out_val = it->second; return true; }
    return false;
}

struct SvgTransform {
    double xmin = 0, xmax = 1, ymin = 0, ymax = 1;
    int width = 800, height = 600;
    double margin = 10.0;
    double scale = 1.0;

    std::pair<double,double> map_xy(const Point_2 &p) const {
        double x = CGAL::to_double(p.x());
        double y = CGAL::to_double(p.y());
        double sx = (xmax > xmin) ? (width - 2.0*margin) / (xmax - xmin) : 1.0;
        double sy = (ymax > ymin) ? (height - 2.0*margin) / (ymax - ymin) : 1.0;
        double s = std::min(sx, sy);
        double xpad = (width - 2.0*margin - s * (xmax - xmin)) * 0.5;
        double ypad = (height - 2.0*margin - s * (ymax - ymin)) * 0.5;
        double X = margin + xpad + (x - xmin) * s;
        double Y = margin + ypad + (ymax - y) * s; // flip y
        return {X, Y};
    }
};

static SvgTransform compute_transform(
    const Arrangement_with_history_2 &arr,
    const std::vector<Point_2> &invalid_points,
    int width = 800,
    int height = 600,
    double margin = 10.0
) {
    double xmin = std::numeric_limits<double>::infinity();
    double ymin = std::numeric_limits<double>::infinity();
    double xmax = -std::numeric_limits<double>::infinity();
    double ymax = -std::numeric_limits<double>::infinity();

    auto update_bbox = [&](const Point_2 &p) {
        double x = CGAL::to_double(p.x());
        double y = CGAL::to_double(p.y());
        xmin = std::min(xmin, x);
        ymin = std::min(ymin, y);
        xmax = std::max(xmax, x);
        ymax = std::max(ymax, y);
    };

    for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
        update_bbox(vit->point());
    }
    for (const auto &p : invalid_points) update_bbox(p);

    if (!std::isfinite(xmin)) {
        xmin = ymin = 0.0;
        xmax = ymax = 1.0;
    }
    if (std::abs(xmax - xmin) < 1e-12) { xmax = xmin + 1.0; }
    if (std::abs(ymax - ymin) < 1e-12) { ymax = ymin + 1.0; }

    SvgTransform T;
    T.xmin = xmin; T.xmax = xmax; T.ymin = ymin; T.ymax = ymax;
    T.width = width; T.height = height; T.margin = margin;
    return T;
}

static bool is_break_vertex(
    const Arrangement_with_history_2::Vertex_const_handle &vh,
    const std::map<Point_2, bool> &point_is_singularity
) {
    bool singular = false;
    auto it = point_is_singularity.find(vh->point());
    if (it != point_is_singularity.end()) singular = it->second;
    int deg = vh->degree();
    return singular || (deg > 2);
}

static bool halfedge_has_single_origin(const Arrangement_with_history_2 &arr,
                                       const Arrangement_with_history_2::Halfedge_const_handle &heh,
                                       Segment_2 &origin_seg) {
    if (arr.number_of_originating_curves(heh) != 1) return false;
    auto ocit = arr.originating_curves_begin(heh);
    origin_seg = *ocit;
    return true;
}

static bool same_origin_segment(const Segment_2 &a, const Segment_2 &b) {
    return same_segment(a, b);
}

static bool is_start_halfedge(
    const Arrangement_with_history_2 &arr,
    const Arrangement_with_history_2::Halfedge_const_handle &heh,
    const Segment_2 &origin_seg,
    const std::map<Point_2, bool> &point_is_singularity
) {
    auto src = heh->source();
    if (is_break_vertex(src, point_is_singularity)) return true;

    auto first = src->incident_halfedges();
    if (first == nullptr) return true;
    auto curr = first;
    do {
        auto h = curr;
        if (h == heh->twin()) {
            ++curr;
            continue;
        }
        if (h->target() == src) {
            Segment_2 s2;
            if (halfedge_has_single_origin(arr, h, s2)) {
                if (same_origin_segment(s2, origin_seg)) {
                    return false;
                }
            }
        }
        ++curr;
    } while (curr != first);

    return true;
}

struct Chain {
    std::vector<Point_2> pts;
    int qi;
    bool convex;
};

std::vector<Chain> get_chains(
    const Arrangement_with_history_2 &arr,
    const std::map<Segment_2, int, Segment2Comparator> &qi,
    const std::map<Segment_2, bool, Segment2Comparator> &segment_is_convex,
    const std::map<Point_2, bool> &point_is_singularity
) {
    std::vector<Chain> chains;
    std::map<Segment_2, int, Segment2Comparator> chain_id;

    for (auto hei = arr.halfedges_begin(); hei != arr.halfedges_end(); ++hei) {
        auto segment = hei->curve();

        if (chain_id.find(segment) != chain_id.end()) {
            continue;
        }

        Chain chain;
        chain.qi = qi.at(segment);

        Segment_2 origin_seg = segment;
        if (arr.number_of_originating_curves(hei) > 0) {
            auto ocit = arr.originating_curves_begin(hei);
            origin_seg = *ocit;
        }

        chain.convex = segment_is_convex.at(origin_seg);

        // traverse the chain
        std::queue<Arrangement_with_history_2::Halfedge_const_handle> queue;
        queue.push(hei);

        std::vector<Segment_2> hes_in_chain;
        std::set<Point_2> pts_in_chain;
        std::map<Point_2, std::vector<Segment_2>> pts_to_segments;
        
        // std::cout << "Traversing chain" << std::endl;

        while (!queue.empty()) {
            auto curr = queue.front();
            queue.pop();

            if (chain_id.find(curr->curve()) != chain_id.end()) {
                continue;
            }

            chain_id[curr->curve()] = chains.size();
            hes_in_chain.push_back(curr->curve());
            pts_in_chain.insert(curr->target()->point());
            pts_in_chain.insert(curr->source()->point());
            pts_to_segments[curr->target()->point()].push_back(curr->curve());
            pts_to_segments[curr->source()->point()].push_back(curr->curve());

            // add the adjacent halfedges to the queue
            std::vector<Arrangement_with_history_2::Vertex_const_handle> adjacent_vertices = {curr->source(), curr->target()};

            for (const auto &v : adjacent_vertices) {
                if (is_break_vertex(v, point_is_singularity)) {
                    continue;
                }

                auto it = v->incident_halfedges();
                do {
                    queue.push(it);
                    ++it;
                } while (it != v->incident_halfedges());
            }
        }

        // std::cout << "Pts in chain: " << pts_in_chain.size() << std::endl;

        Point_2 start_pt;
        for (const auto &p : pts_in_chain) {
            if (pts_to_segments[p].size() == 1) {
                start_pt = p;
                break;
            }
        }
        // std::cout << "Start point 1: " << start_pt << std::endl;

        if (start_pt == Point_2()) {
            // if no start point is found, then any point is okay for start_pt
            // but we need to find one that actually has segments
            for (const auto &p : pts_in_chain) {
                if (pts_to_segments[p].size() > 0) {
                    start_pt = p;
                    break;
                }
            }
        }

        // std::cout << "Start point 2: " << start_pt << std::endl;

        // Check if we found a valid start point with segments
        if (start_pt == Point_2() || pts_to_segments[start_pt].empty()) {
            std::cout << "Warning: No valid start point found for chain, skipping..." << std::endl;
            continue;
        }

        std::map<Point_2, bool> pts_visited;

        chain.pts.push_back(start_pt);
        pts_visited[start_pt] = true;

        Segment_2 current_seg = pts_to_segments[start_pt][0];
        Point_2 current_pt = current_seg.target() == start_pt ? current_seg.source() : current_seg.target();

        while (true) {
            chain.pts.push_back(current_pt);
            pts_visited[current_pt] = true;

            auto segments = pts_to_segments[current_pt];
            if (segments.size() == 1) {
                // std::cout << "Breaking at point due to single segment: " << current_pt << std::endl;
                break;
            }

            current_seg = same_segment(current_seg, segments[0]) ? segments[1] : segments[0];
            current_pt = current_seg.target() == current_pt ? current_seg.source() : current_seg.target();

            if (pts_visited.find(current_pt) != pts_visited.end()) {
                // std::cout << "Breaking at point due to loop: " << current_pt << std::endl;
                // since it's a loop, add the last points
                chain.pts.push_back(current_pt);
                break;
            }
        }

        // std::cout << "Chain points: " << chain.pts.size() << std::endl;

        if (chain.pts.size() != pts_in_chain.size() && chain.pts.size() != pts_in_chain.size() + 1) {
            std::cout << "Chain points: " << chain.pts.size() << " Pts in chain: " << pts_in_chain.size() << std::endl;
            throw std::runtime_error("Chain points and pts in chain size mismatch");
        }

        chains.push_back(chain);
    }

    return chains;
}

std::string planar_map_to_svg(
    const Arrangement_with_history_2 &arr,
    const std::map<Segment_2, int, Segment2Comparator> &qi,
    const std::map<Segment_2, bool, Segment2Comparator> &segment_is_convex,
    const std::map<Point_2, bool> &point_is_singularity,
    const std::vector<Point_2> &invalid_points
) {
    SvgTransform T = compute_transform(arr, invalid_points, 800, 600, 10.0);

    auto chains = get_chains(arr, qi, segment_is_convex, point_is_singularity);

    std::ostringstream ss;
    ss << "<svg xmlns=\"http://www.w3.org/2000/svg\""
       << " width=\"" << T.width << "\" height=\"" << T.height << "\""
       << " viewBox=\"0 0 " << T.width << " " << T.height << "\">" << "\n";

    ss << "  <g id=\"arrangement\">\n";
    ss << std::fixed << std::setprecision(3);

    for (const auto &ch : chains) {
        if (ch.pts.size() < 2) continue;

        std::string color = color_for_qi(ch.qi);
        ss << "    <path d=\"";
        for (size_t i = 0; i < ch.pts.size(); ++i) {
            auto xy = T.map_xy(ch.pts[i]);
            double x = xy.first, y = xy.second;
            if (i == 0) {
                ss << "M " << x << " " << y << " ";
            } else {
                ss << "L " << x << " " << y << " ";
            }
        }
        ss << "\" fill=\"none\""
           << " stroke=\"" << color << "\""
           << " stroke-width=\"" << 1.0 << "\"";
        if (!ch.convex) {
            ss << " stroke-dasharray=\"1 2\" stroke-linecap=\"round\"";
        }
        ss << "/>\n";
    }
    ss << "  </g>\n";

    ss << "  <g id=\"singularities\">\n";
    for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
        auto pit = point_is_singularity.find(vit->point());
        if (pit != point_is_singularity.end() && pit->second) {
            auto xy = T.map_xy(vit->point());
            double x = xy.first, y = xy.second;
            ss << "    <circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"1\" fill=\"black\" stroke=\"none\"/>\n";
        }
    }
    ss << "  </g>\n";

    ss << "  <g id=\"invalid_points\">\n";
    for (const auto &p : invalid_points) {
        auto xy = T.map_xy(p);
        double x = xy.first, y = xy.second;
        ss << "    <circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"1.5\" fill=\"none\" stroke=\"red\"/>\n";
    }
    ss << "  </g>\n";

    ss << "</svg>\n";
    return ss.str();
}
}
