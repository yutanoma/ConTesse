#include "tessellate_face.h"

#include "igl/edge_topology.h"
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_with_history_2.h>
#include <vector>
#include <map>

namespace utils {

std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXi,
           Eigen::MatrixXi, Eigen::MatrixXi, std::map<Point_2, int>,
           std::map<Segment_2, int, Segment2Comparator>>
tessellate_face(const Arrangement_with_history_2 &arr,
                const Face_const_handle &face)
{
  // Type definitions for constrained Delaunay triangulation
  typedef CGAL::Triangulation_vertex_base_2<Kernel> Vb;
  typedef CGAL::Constrained_triangulation_face_base_2<Kernel> Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
  typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, Tds> CDT;
  typedef CDT::Point CDT_Point;
  typedef CDT::Vertex_handle Vertex_handle;

  // Extract face boundary points and constraints
  std::vector<CDT_Point> points;
  std::vector<std::pair<int, int>> constraints;
  std::map<Point_2, int> point_to_vid;
  std::map<Segment_2, int, Segment2Comparator> segment_to_eid;
  
  int vertex_id = 0;
  
  // Helper function to add a point and return its ID
  auto add_point = [&](const Point_2& p) -> int {
    auto it = point_to_vid.find(p);
    if (it != point_to_vid.end()) {
      return it->second;
    }
    
    int id = vertex_id++;
    point_to_vid[p] = id;
    points.push_back(CDT_Point(CGAL::to_double(p.x()), CGAL::to_double(p.y())));
    return id;
  };
  
  // Helper function to add a constraint edge
  auto add_constraint = [&](int v1, int v2) {
    constraints.push_back({v1, v2});
  };

  // Process outer boundary
  if (face->has_outer_ccb()) {
    auto outer_ccb = face->outer_ccb();
    if (outer_ccb != CGAL::Arrangement_2<Traits_2>::Ccb_halfedge_circulator()) {
      std::vector<int> outer_vertices;
      auto he = outer_ccb;
      do {
        Point_2 source = he->source()->point();
        int vid = add_point(source);
        outer_vertices.push_back(vid);
        ++he;
      } while (he != outer_ccb);
      
      // Add constraints for outer boundary
      for (size_t i = 0; i < outer_vertices.size(); ++i) {
        int next = (i + 1) % outer_vertices.size();
        add_constraint(outer_vertices[i], outer_vertices[next]);
      }
    }
  }
  
  // Process holes (inner boundaries)
  for (auto hole_it = face->holes_begin(); hole_it != face->holes_end(); ++hole_it) {
    auto hole_ccb = *hole_it;
    std::vector<int> hole_vertices;
    auto he = hole_ccb;
    do {
      Point_2 source = he->source()->point();
      int vid = add_point(source);
      hole_vertices.push_back(vid);
      ++he;
    } while (he != hole_ccb);
    
    // Add constraints for hole boundary
    for (size_t i = 0; i < hole_vertices.size(); ++i) {
      int next = (i + 1) % hole_vertices.size();
      add_constraint(hole_vertices[i], hole_vertices[next]);
    }
  }

  // Create constrained Delaunay triangulation
  CDT cdt;
  
  // Insert points
  for (const auto& point : points) {
    cdt.insert(point);
  }
  
  // Insert constraints
  for (const auto& constraint : constraints) {
    CDT_Point p1 = points[constraint.first];
    CDT_Point p2 = points[constraint.second];
    cdt.insert_constraint(p1, p2);
  }

  // Extract vertices and faces from triangulation
  std::map<Vertex_handle, int> vertex_to_id;
  int current_id = 0;
  
  // Build vertex mapping
  for (auto vit = cdt.vertices_begin(); vit != cdt.vertices_end(); ++vit) {
    vertex_to_id[vit] = current_id++;
  }
  
  // Create vertex matrix
  Eigen::MatrixXd V(current_id, 2);
  for (auto vit = cdt.vertices_begin(); vit != cdt.vertices_end(); ++vit) {
    int id = vertex_to_id[vit];
    V(id, 0) = CGAL::to_double(vit->point().x());
    V(id, 1) = CGAL::to_double(vit->point().y());
  }
  
  // Create face matrix - include all faces from the triangulation
  // The constrained Delaunay triangulation will only create faces within the constraints
  std::vector<Eigen::Vector3i> faces;
  for (auto fit = cdt.faces_begin(); fit != cdt.faces_end(); ++fit) {
    // Skip the infinite face
    if (cdt.is_infinite(fit)) {
      continue;
    }
    
    int v0 = vertex_to_id[fit->vertex(0)];
    int v1 = vertex_to_id[fit->vertex(1)];
    int v2 = vertex_to_id[fit->vertex(2)];
    faces.push_back(Eigen::Vector3i(v0, v1, v2));
  }
  
  Eigen::MatrixXi F(faces.size(), 3);
  for (size_t i = 0; i < faces.size(); ++i) {
    F.row(i) = faces[i];
  }

  // Compute edge topology using igl
  Eigen::MatrixXi EV, FE, EF;
  igl::edge_topology(V, F, EV, FE, EF);

  // Now map segments to actual edge IDs from EV, FE, EF
  // Create a mapping from vertex pairs to edge IDs
  std::map<std::pair<int, int>, int> vertex_pair_to_edge_id;
  for (int e = 0; e < EV.rows(); ++e) {
    int v1 = EV(e, 0);
    int v2 = EV(e, 1);
    // Store both directions for easy lookup
    vertex_pair_to_edge_id[{v1, v2}] = e;
    vertex_pair_to_edge_id[{v2, v1}] = e;
  }
  
  // Map original segments to edge IDs
  for (const auto& constraint : constraints) {
    int v1 = constraint.first;
    int v2 = constraint.second;
    
    // Find the corresponding edge ID
    auto it = vertex_pair_to_edge_id.find({v1, v2});
    if (it != vertex_pair_to_edge_id.end()) {
      int edge_id = it->second;
      
      // Create the segment from the original points
      Point_2 p1 = points[v1];
      Point_2 p2 = points[v2];
      Segment_2 seg(p1, p2);
      segment_to_eid[seg] = edge_id;
    }
  }

  return {V, F, EV, FE, EF, point_to_vid, segment_to_eid};
}
}
