#include "planar_map_to_svg.h"
#include "planar_map.h"

#include <sstream>
#include <iomanip>
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
};

static SvgTransform compute_transform(
    const Arrangement_with_history_2 &arr,
    const std::vector<Point_2> &invalid_points
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
        xmin = ymin = 0;
        xmax = ymax = 800;
    }
    if (std::abs(xmax - xmin) < 1e-12) { xmax = xmin + 1.0; }
    if (std::abs(ymax - ymin) < 1e-12) { ymax = ymin + 1.0; }

    SvgTransform T;
    T.xmin = xmin;
    T.xmax = xmax;
    T.ymin = ymin;
    T.ymax = ymax;

    T.width = xmax;
    T.height = ymax;

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

struct Chain {
    std::vector<Point_2> pts;
    int qi;
    bool convex;
    bool all_cut;
};

std::vector<Chain> get_chains(
    const Arrangement_with_history_2 &arr,
    const std::map<Segment_2, int, Segment2Comparator> &qi,
    const std::map<Segment_2, bool, Segment2Comparator> &segment_is_convex,
    const std::map<Segment_2, bool, Segment2Comparator> &segment_is_cut,
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
        chain.convex = segment_is_convex.at(segment);
        {
            auto it_cut_init = segment_is_cut.find(segment);
            bool is_cut_init = (it_cut_init != segment_is_cut.end()) ? it_cut_init->second : false;
            chain.all_cut = is_cut_init; // chain class is defined by first segment's cut state
        }

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
                    // Only expand along segments that match the chain's cut classification
                    if (chain_id.find(it->curve()) == chain_id.end()) {
                        auto it_cut_n = segment_is_cut.find(it->curve());
                        bool is_cut_n = (it_cut_n != segment_is_cut.end()) ? it_cut_n->second : false;
                        if (is_cut_n == chain.all_cut) {
                            queue.push(it);
                        }
                    }
                    ++it;
                } while (it != v->incident_halfedges());
            }
        }

        // std::cout << "Pts in chain: " << pts_in_chain.size() << std::endl;

        // std::cout << "default point: " << Point_2() << std::endl;

        Point_2 start_pt;
        bool found_start_pt = false;
        for (const auto &p : pts_in_chain) {
            if (pts_to_segments[p].size() == 1) {
                start_pt = p;
                found_start_pt = true;
                break;
            }
        }

        if (!found_start_pt) {
            // if no start point is found, then any point is okay for start_pt
            // but we need to find one that actually has segments
            for (const auto &p : pts_in_chain) {
                if (pts_to_segments[p].size() > 0) {
                    start_pt = p;
                    found_start_pt = true;
                    break;
                }
            }
        }

        // Check if we found a valid start point with segments
        if (!found_start_pt || pts_to_segments[start_pt].empty()) {
            std::cout << "Warning: No valid start point found for chain, skipping..." << std::endl;
            continue;
        }

        std::map<Point_2, bool> pts_visited;

        chain.pts.push_back(start_pt);
        pts_visited[start_pt] = true;

        Segment_2 current_seg = pts_to_segments[start_pt][0];
        Point_2 current_pt = is_identical_point(current_seg.target(), start_pt) ? current_seg.source() : current_seg.target();

        while (true) {
            chain.pts.push_back(current_pt);
            pts_visited[current_pt] = true;

            auto segments = pts_to_segments[current_pt];
            if (segments.size() == 1) {
                // std::cout << "Breaking at point due to single segment: " << current_pt << std::endl;
                break;
            }

            current_seg = is_identical_segment(current_seg, segments[0]) ? segments[1] : segments[0];
            current_pt = is_identical_point(current_seg.target(), current_pt) ? current_seg.source() : current_seg.target();

            if (pts_visited.find(current_pt) != pts_visited.end()) {
                // std::cout << "Breaking at point due to loop: " << current_pt << std::endl;
                // since it's a loop, add the last points
                chain.pts.push_back(current_pt);
                break;
            }
        }

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
    const std::map<Segment_2, bool, Segment2Comparator> &segment_is_cut,
    const std::map<Point_2, bool> &point_is_singularity,
    const std::vector<Point_2> &invalid_points
) {
    SvgTransform T = compute_transform(arr, invalid_points);

    auto chains = get_chains(arr, qi, segment_is_convex, segment_is_cut, point_is_singularity);

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
            double x = CGAL::to_double(ch.pts[i].x());
            double y = CGAL::to_double(ch.pts[i].y());
            if (i == 0) {
                ss << "M " << x << " " << y << " ";
            } else {
                ss << "L " << x << " " << y << " ";
            }
        }
        ss << "\" fill=\"none\""
           << " stroke=\"" << color << "\""
           << " stroke-width=\"" << 1.0 << "\"";
        if (ch.all_cut) {
            ss << " stroke-dasharray=\"4 6\" stroke-linecap=\"round\"";
        } else if (!ch.convex) {
            ss << " stroke-dasharray=\"1 2\" stroke-linecap=\"round\"";
        }
        ss << "/>\n";
    }
    ss << "  </g>\n";

    ss << "  <g id=\"singularities\">\n";
    for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
        auto pit = point_is_singularity.find(vit->point());
        if (pit != point_is_singularity.end() && pit->second) {
            double x = CGAL::to_double(vit->point().x());
            double y = CGAL::to_double(vit->point().y());
            ss << "    <circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"1\" fill=\"black\" stroke=\"none\"/>\n";
        }
    }
    ss << "  </g>\n";

    ss << "  <g id=\"invalid_points\">\n";
    for (const auto &p : invalid_points) {
        double x = CGAL::to_double(p.x());
        double y = CGAL::to_double(p.y());
        ss << "    <circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"1.5\" fill=\"none\" stroke=\"red\"/>\n";
    }
    ss << "  </g>\n";

    ss << "</svg>\n";
    return ss.str();
}
}
