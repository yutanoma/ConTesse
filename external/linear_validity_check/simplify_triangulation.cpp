#include "simplify_triangulation.h"

#include <Eigen/Core>
#include <igl/triangle_triangle_adjacency.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <vector>
#include <unordered_set>
#include <queue>

namespace utils {
namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;

using d2       = bg::model::d2::point_xy<double>;
using Polygon  = bg::model::polygon<d2, /*Clockwise=*/true, /*Closed=*/true>;
using Box      = bg::model::box<d2>;

// Helper: make closed triangle polygon
static inline Polygon make_triangle(const Eigen::MatrixXd& V, int a, int b, int c) {
    Polygon tri;
    auto& out = tri.outer();
    out.reserve(4);
    out.emplace_back(V(a,0), V(a,1));
    out.emplace_back(V(b,0), V(b,1));
    out.emplace_back(V(c,0), V(c,1));
    out.emplace_back(V(a,0), V(a,1));
    bg::correct(tri);
    return tri;
}

// Compute strict positive-area intersection test
static inline bool polygons_overlap_with_area(const Polygon& A, const Polygon& B, double eps = 0.0) {
    std::vector<Polygon> inter;
    bg::intersection(A, B, inter);
    double area = 0.0;
    for (const auto& p : inter) area += bg::area(p);
    return area > eps;
}

// Main function
std::vector<std::vector<int>> overlapping_faces_2d_bvh(const Eigen::MatrixXd &V,
                                                       const Eigen::MatrixXi &F,
                                                       const Eigen::MatrixXi &TT,
                          bool strict_positive_area = false,
                          bool indices_are_one_based = false)
{
    using Value = std::pair<Box,int>;
    const int m = static_cast<int>(F.rows());

    // Store adjacency sets for quick lookup
    std::vector<std::unordered_set<int>> neighbors(m);
    for (int i = 0; i < m; ++i) {
        for (int k = 0; k < 3; ++k) {
            int nb = TT(i, k);
            if (nb >= 0) neighbors[i].insert(nb);
        }
    }

    // Prepare polygons + bounding boxes
    auto idx = [&](int x)->int { return indices_are_one_based ? (x-1) : x; };
    std::vector<Polygon> tris(m);
    std::vector<Value> items;
    items.reserve(m);
    for (int i = 0; i < m; ++i) {
        int a = idx(F(i,0)), b = idx(F(i,1)), c = idx(F(i,2));
        tris[i] = make_triangle(V, a, b, c);
        Box box; bg::envelope(tris[i], box);
        items.emplace_back(box, i);
    }

    // Build BVH
    bgi::rtree<Value, bgi::rstar<16>> rtree(items.begin(), items.end());

    // Query
    std::vector<std::vector<int>> overlaps(m);
    std::vector<Value> candidates;

    for (int i = 0; i < m; ++i) {
        const Polygon& Pi = tris[i];
        candidates.clear();
        rtree.query(bgi::intersects(items[i].first), std::back_inserter(candidates));

        for (const auto& [boxj, j] : candidates) {
            if (j == i) continue;
            if (neighbors[i].count(j)) continue; // skip incident faces

            const Polygon& Pj = tris[j];
            bool hit = false;
            if (!strict_positive_area) {
                if (bg::intersects(Pi, Pj)) hit = true;
            } else {
                if (polygons_overlap_with_area(Pi, Pj)) hit = true;
            }
            if (hit) overlaps[i].push_back(j);
        }
    }

    return overlaps;
}

// Group faces so that within each group, no two faces overlap (per `overlaps`).
// Expansion is along incident faces given by TT (triangle-triangle adjacency).
// Returns group id per face (0,1,2,...) and -1 for meshes with zero faces.
Eigen::VectorXi
group_triangulation(
    const Eigen::MatrixXd& V,                        // n x 2 (not used here; kept for signature parity)
    const Eigen::MatrixXi& F,                        // m x 3
    const Eigen::MatrixXi& TT,                       // m x 3, -1 for no neighbor
    const std::vector<std::vector<int>>& overlaps    // size m; overlaps[i] are face ids that overlap face i
) {
    const int m = static_cast<int>(F.rows());
    Eigen::VectorXi group_of = Eigen::VectorXi::Constant(m, -1);
    if (m == 0) return group_of;

    // Optional: pre-build an "overlaps bit" structure per face for faster conflict checks.
    // For simplicity and memory efficiency, we’ll check against group members’ marker.
    // We’ll maintain a boolean array `in_current_group` using a generation counter to avoid O(m) clears.
    std::vector<int> in_group_stamp(m, -1);
    int current_stamp = 0;

    int gid = 0;

    // Iterate over all faces as potential seeds
    for (int seed = 0; seed < m; ++seed) {
        if (group_of(seed) != -1) continue; // already placed in some group

        // Start a new group from this seed
        std::queue<int> q;          // work queue within this group
        std::vector<int> members;   // list of faces assigned to this group (for final stamping/checks)

        // Start new generation of in-group markers
        ++current_stamp;

        // Accept the seed into the group immediately
        q.push(seed);
        group_of(seed) = gid;
        in_group_stamp[seed] = current_stamp;
        members.push_back(seed);

        // BFS along incident faces
        while (!q.empty()) {
            int f = q.front(); q.pop();

            // Explore 3 incident neighbors of face f via TT
            for (int k = 0; k < 3; ++k) {
                int nb = TT(f, k);
                if (nb < 0) continue;
                if (group_of(nb) != -1) continue; // already assigned to some group

                // Check conflict: does `nb` overlap with any face already in this group?
                bool conflict = false;
                for (int o : overlaps[nb]) {
                    if (o >= 0 && o < m && in_group_stamp[o] == current_stamp) {
                        conflict = true;
                        break;
                    }
                }
                if (conflict) continue;

                // Safe to add `nb` into this group
                group_of(nb) = gid;
                in_group_stamp[nb] = current_stamp;
                members.push_back(nb);
                q.push(nb);
            }
        }

        // Done filling this group; advance group id
        ++gid;
    }

    return group_of;
}

// Packs (group, vertex) into a 64-bit key (requires group, vertex >= 0)
static inline uint64_t pack_key(int group, int v) {
  return (static_cast<uint64_t>(static_cast<uint32_t>(group)) << 32) |
         (static_cast<uint64_t>(static_cast<uint32_t>(v)));
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXi>
unify_simplified_meshes(
    const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
    const Eigen::VectorXi &group_id
) {
    const int n = static_cast<int>(V.rows());
    const int m = static_cast<int>(F.rows());
    assert(F.cols() == 3 && "F must be m x 3");
    assert(V.cols() == 2 && "V must be n x 2 (2D)");
    assert(group_id.size() == m && "group_id must be length m");

    // Map (group, original vertex) -> new vertex index
    std::unordered_map<uint64_t, int> gv2new;
    gv2new.reserve(static_cast<size_t>(m) * 3);

    // For output
    std::vector<int> Vnew_to_orig;
    Eigen::MatrixXi Fnew(m, 3);

    auto ensure_new = [&](int g, int v)->int {
        assert(g >= 0 && "group_id must be non-negative");
        assert(v >= 0 && v < n && "vertex index out of range");
        const uint64_t key = pack_key(g, v);
        auto it = gv2new.find(key);
        if (it != gv2new.end()) return it->second;

        const int new_id = static_cast<int>(Vnew_to_orig.size());
        gv2new.emplace(key, new_id);
        Vnew_to_orig.push_back(v);
        return new_id;
    };

    // Build F' by assigning per-(group,origV) vertices
    for (int f = 0; f < m; ++f) {
        const int g = group_id(f);
        assert(g >= 0 && "group_id(f) must be >= 0");
        const int v0 = F(f,0);
        const int v1 = F(f,1);
        const int v2 = F(f,2);
        const int nv0 = ensure_new(g, v0);
        const int nv1 = ensure_new(g, v1);
        const int nv2 = ensure_new(g, v2);
        Fnew(f,0) = nv0;
        Fnew(f,1) = nv1;
        Fnew(f,2) = nv2;
    }

    // Finalize V'
    Eigen::MatrixXd Vnew(static_cast<int>(Vnew_to_orig.size()), 2);
    for (int i = 0; i < static_cast<int>(Vnew_to_orig.size()); ++i) {
        Vnew(i,0) = V(Vnew_to_orig[i], 0);
        Vnew(i,1) = V(Vnew_to_orig[i], 1);
    }

    Eigen::VectorXi Vnew_to_orig_vec =
        Eigen::Map<Eigen::VectorXi>(Vnew_to_orig.data(), Vnew_to_orig.size());

    return {Vnew, Fnew, Vnew_to_orig_vec};
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXi>
simplify_triangulation(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F) {
  // 1. get the overlapping faces
  Eigen::MatrixXi TT, TTi;
  igl::triangle_triangle_adjacency(F, TT, TTi);
  auto overlapping_faces = overlapping_faces_2d_bvh(V, F, TT);

  // 2. based on the overlapping faces, simplify the triangulation
  auto face_groups = group_triangulation(V, F, TT, overlapping_faces);

  // 3. unify it as a single triangle mesh
  auto [V_simp, F_simp, V_simp_to_out] = unify_simplified_meshes(V, F, face_groups);
  
  return {V_simp, F_simp, V_simp_to_out};
}
}
