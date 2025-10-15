#pragma once

#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_with_history_2.h>
#include <Eigen/Core>
#include <vector>
#include <map>
#include <unordered_map>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel> Traits_2;
typedef Traits_2::Point_2 Point_2;
typedef Traits_2::X_monotone_curve_2 Segment_2;
typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;
typedef CGAL::Arrangement_with_history_2<Traits_2> Arrangement_with_history_2;
typedef CGAL::Arrangement_with_history_2<Traits_2>::Face_const_handle Face_const_handle;
typedef CGAL::Arrangement_with_history_2<Traits_2>::Halfedge_const_handle Halfedge_const_handle;

// Custom hash and equality for CGAL keys
struct Point2Hash {
    size_t operator()(const Point_2& p) const noexcept {
        size_t h1 = std::hash<double>()(CGAL::to_double(p.x()));
        size_t h2 = std::hash<double>()(CGAL::to_double(p.y()));
        return h1 ^ (h2 << 1);
    }
};
struct Point2Equal {
    bool operator()(const Point_2& a, const Point_2& b) const noexcept {
        return a == b;
    }
};
struct Segment2Hash {
    size_t operator()(const Segment_2& s) const noexcept {
        size_t h1 = Point2Hash{}(s.source());
        size_t h2 = Point2Hash{}(s.target());
        return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
    }
};
struct Segment2Equal {
    bool operator()(const Segment_2& a, const Segment_2& b) const noexcept {
        return a.source() == b.source() && a.target() == b.target();
    }
};

template <typename T>
using Point2List = std::unordered_map<Point_2, T, Point2Hash, Point2Equal>;

template <typename T>
using Segment2List = std::unordered_map<Segment_2, T, Segment2Hash, Segment2Equal>;

namespace utils {

// Custom comparator for Segment_2 objects
struct Segment2Comparator {
    bool operator()(const Segment_2& a, const Segment_2& b) const {
        auto a_pair = std::make_pair(a.source(), a.target());
        auto b_pair = std::make_pair(b.source(), b.target());
        return a_pair < b_pair;
    }
};

bool is_identical_segment(const Segment_2 &s1, const Segment_2 &s2);

std::tuple<
    Arrangement_with_history_2,
    std::map<Point_2, int>,
    std::map<Segment_2, int, Segment2Comparator>
> planar_map(
    const Eigen::MatrixXd &V_2d,
    const Eigen::MatrixXi &E_connectivity
);
}
