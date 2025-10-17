#include "linear_wso_check.h"

#include "common.h"
#include "triangulate_wso_patches.h"

#include "linear_validity_check/fast_validity_check.h"

#include <fstream>
#include <filesystem>

void linear_wso_check(Mesh &mesh, Camera const &camera,
                      std::map<size_t, Mesh> &patch_triangulations,
                      bool do_refine,
                      std::unordered_set<int> selected_patches,
                      bool ff_only
) {
  auto is_cusp = mesh.get_vertex_property<bool>("v:cusp");
  auto cusp_facing = mesh.get_vertex_property<FacingType>("v:cusp_facing");

  auto cut = mesh.get_edge_property<int>("e:disk_cut");
  contess_assert_msg(cut, "Cuts are needed for chaining cut graph.");

  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  Vector2i viewport({camera.vpWidth(), camera.vpHeight()});

  // Make sure the mesh has patches extracted
  contess_assert_msg(!mesh.get_seams_lengths_const().empty(),
                     "Can't triangulate when no patch is extracted.");
  // Make sure the rotation indices are computed
  contess_assert_msg(is_cusp && cusp_facing,
                     "Can't triangulate when no cusp is extracted.");

  patch_triangulations.clear();

  auto vertex_to_idx = mesh.add_vertex_property<int>("v:wso_idx", -1);
  // Initialize all values to -1
  vertex_to_idx.vector().assign(vertex_to_idx.vector().size(), -1);

  auto edge_to_idx = mesh.add_edge_property<int>("e:wso_idx", -1);
  // Initialize all values to -1
  edge_to_idx.vector().assign(edge_to_idx.vector().size(), -1);

  // Remove existing debug inputs directory if it exists
  {
    const std::string root_debug_dir = "./debug_inputs";
    if (std::filesystem::exists(root_debug_dir)) {
      std::filesystem::remove_all(root_debug_dir);
    }
  }

  for (auto const &patch_contour : mesh.get_const_patch_chains()) {
    int patch_id = patch_contour.first;

    if (selected_patches.empty() || selected_patches.count(patch_id))
      logger().info("WSO on {}.", patch_id);
    if (!selected_patches.empty() && !selected_patches.count(patch_id)) {
      continue;
    }

    FacingType facing = mesh.get_patch_facing(patch_id);

    if (ff_only && facing != FacingType::FRONT)
      continue;

    // it's a cut if the chain traverses the edge twice
    // otherwise, it's a boundary
    // mark the orientation and the cut
    for (auto const &he : *patch_contour.second) {
      logger().info("Halfedge: {} -> {}", mesh.from_vertex(he).idx(), mesh.to_vertex(he).idx());

      // first, mark all the vertices in the chain
      vertex_to_idx[mesh.from_vertex(he)] = -1;
      vertex_to_idx[mesh.to_vertex(he)] = -1;
    }

    int idx = 0;
    std::vector<Vertex> vertices;
    for (auto const &he : *patch_contour.second) {
      if (vertex_to_idx[mesh.from_vertex(he)] == -1) {
        vertex_to_idx[mesh.from_vertex(he)] = idx;
        vertices.push_back(mesh.from_vertex(he));
        idx++;
      }
      if (vertex_to_idx[mesh.to_vertex(he)] == -1) {
        vertex_to_idx[mesh.to_vertex(he)] = idx;
        vertices.push_back(mesh.to_vertex(he));
        idx++;
      }
    }

    int e_idx = 0;
    std::vector<Edge> edges;
    std::vector<int> edge_signs;
    std::vector<bool> is_edge_cut;
    for (auto const &he : *patch_contour.second) {
      if (edge_to_idx[mesh.edge(he)] == -1) {
        edge_to_idx[mesh.edge(he)] = e_idx;
        edges.push_back(mesh.edge(he));

        if (cut[mesh.edge(he)] < 0) {
          is_edge_cut.push_back(true);
          edge_signs.push_back(1);

          logger().info("Edge {} is a cut", mesh.edge(he).idx());
        } else {
          is_edge_cut.push_back(false);

          // this is a boundary edge
          // the sign is determined by the orientation of the chain
          auto edge = mesh.edge(he);
          auto v_from = mesh.vertex(edge, 0);
          auto v_to = mesh.vertex(edge, 1);

          if (v_from == mesh.from_vertex(he) && v_to == mesh.to_vertex(he)) {
            edge_signs.push_back(1);
          } else if (v_from == mesh.to_vertex(he) && v_to == mesh.from_vertex(he)) {
            edge_signs.push_back(-1);
          } else {
            logger().error("Boundary edge is not consistent with the chain orientation.");
            continue;
          }
        }

        e_idx++;
      } else {
        // it's the second time to see this
        // this means that the edge is cut
        // is_edge_cut[edge_to_idx[mesh.edge(he)]] = true;
        logger().info("Edge {}, index is {} is a cut", mesh.edge(he).idx(),
                      edge_to_idx[mesh.edge(he)]);
        is_edge_cut[edge_to_idx[mesh.edge(he)]] = true;
      }
    }

    // if this is a back facing patch, reverse the signs of the edges
    bool is_back_facing = facing == FacingType::BACK;

    Matrix3Xf _vertices_3d(3, vertices.size());
    for (size_t i = 0; i < vertices.size(); i++) {
      _vertices_3d.col(i) = vpositions[vertices[i]];
    }

    std::vector<Vector2f> polyline;

    projectToViewport(_vertices_3d, polyline, camera.projectionMatrix(),
                      camera.viewMatrix(), viewport);

    Eigen::MatrixXd vertices_2d(vertices.size(), 2);
    for (size_t i = 0; i < vertices.size(); i++) {
      vertices_2d(i, 0) = static_cast<double>(polyline[i].x());
      vertices_2d(i, 1) = static_cast<double>(polyline[i].y());
    }

    Eigen::MatrixXi edges_connectivity(edges.size(), 2);
    for (size_t i = 0; i < edges.size(); i++) {
      edges_connectivity(i, 0) = vertex_to_idx[mesh.vertex(edges[i], 0)];
      edges_connectivity(i, 1) = vertex_to_idx[mesh.vertex(edges[i], 1)];
    }

    Eigen::MatrixXd vertices_3d(vertices.size(), 3);
    for (size_t i = 0; i < vertices.size(); i++) {
      for (size_t j = 0; j < 3; j++) {
        vertices_3d(i, j) = static_cast<double>(_vertices_3d(j, i));
      }
    }

    Eigen::VectorXi edges_sign(edges.size());
    for (size_t i = 0; i < edges.size(); i++) {
      edges_sign(i) = edge_signs[i];
    }

    Eigen::VectorXi edges_is_cut(edges.size());
    for (size_t i = 0; i < edges.size(); i++) {
      edges_is_cut(i) = is_edge_cut[i];
    }

    Eigen::VectorXi vertex_is_cusp(vertices.size());
    for (size_t i = 0; i < vertices.size(); i++) {
      if (is_cusp[vertices[i]]) {
        vertex_is_cusp(i) = 1;
        logger().info("Vertex {} is a cusp", vertices[i].idx());
        logger().info("Facing type: {}", cusp_facing[vertices[i]]);
      } else {
        vertex_is_cusp(i) = 0;
      }
    }

    Eigen::VectorXi vertex_is_singularity(vertices.size());
    for (size_t i = 0; i < vertices.size(); i++) {
      if (is_cusp[vertices[i]] && cusp_facing[vertices[i]] != FacingType::UNDEFINED && facing != cusp_facing[vertices[i]]) {
        vertex_is_singularity(i) = 1;
      } else {
        vertex_is_singularity(i) = 0;
      }
    }

    auto camera_pos = static_cast<Eigen::Vector3d>(camera.position());

    // Export inputs for debugging
    std::string debug_dir = "./debug_inputs/" + std::to_string(patch_id);
    std::filesystem::create_directories(debug_dir);
    
    std::ofstream vertices_3d_file(debug_dir + "/vertices_3d.txt", std::ios::out | std::ios::trunc);
    vertices_3d_file << vertices_3d << std::endl;
    vertices_3d_file.close();
    
    std::ofstream vertices_2d_file(debug_dir + "/vertices_2d.txt", std::ios::out | std::ios::trunc);
    vertices_2d_file << vertices_2d << std::endl;
    vertices_2d_file.close();
    
    std::ofstream edges_connectivity_file(debug_dir + "/edges_connectivity.txt", std::ios::out | std::ios::trunc);
    edges_connectivity_file << edges_connectivity << std::endl;
    edges_connectivity_file.close();
    
    std::ofstream edges_sign_file(debug_dir + "/edges_sign.txt", std::ios::out | std::ios::trunc);
    edges_sign_file << edges_sign << std::endl;
    edges_sign_file.close();
    
    std::ofstream edges_is_cut_file(debug_dir + "/edges_is_cut.txt", std::ios::out | std::ios::trunc);
    edges_is_cut_file << edges_is_cut << std::endl;
    edges_is_cut_file.close();
    
    std::ofstream vertex_is_cusp_file(debug_dir + "/vertex_is_cusp.txt", std::ios::out | std::ios::trunc);
    vertex_is_cusp_file << vertex_is_cusp << std::endl;
    vertex_is_cusp_file.close();
    
    std::ofstream vertex_is_singularity_file(debug_dir + "/vertex_is_singularity.txt", std::ios::out | std::ios::trunc);
    vertex_is_singularity_file << vertex_is_singularity << std::endl;
    vertex_is_singularity_file.close();
    
    std::ofstream camera_pos_file(debug_dir + "/camera_pos.txt", std::ios::out | std::ios::trunc);
    camera_pos_file << camera_pos.transpose() << std::endl;
    camera_pos_file.close();
    
    std::ofstream is_back_facing_file(debug_dir + "/is_back_facing.txt", std::ios::out | std::ios::trunc);
    is_back_facing_file << is_back_facing << std::endl;
    is_back_facing_file.close();
    
    logger().info("Debug inputs exported to {} directory", debug_dir);

    logger().info("fast_validity_check start");

    // auto [is_wso_succeeded, V_out, F_out, V_JI] = utils::fast_validity_check(
    //     vertices_3d, vertices_2d, edges_connectivity, edges_sign, edges_is_cut,
    //     vertex_is_cusp, vertex_is_singularity, camera_pos, is_back_facing);

    // std::vector<size_t> vid_original(V_JI.size());
    // for (int i = 0; i < V_JI.size(); i++) {
    //   vid_original[i] = vertices[V_JI(i)].idx();
    // }

    // patch_triangulations[patch_id] = Mesh();
    // to_2d_mesh(mesh, V_out, F_out, vid_original,
    //            patch_triangulations[patch_id]);

    // auto patch_patchID =
    //     patch_triangulations[patch_id].face_property<int>("f:patchID");
    // patch_patchID.vector().assign(patch_patchID.vector().size(), patch_id);

    // if (!is_wso_succeeded) {
    //   logger().error("Linear WSO check failed on patch {}.", patch_id);
    //   continue;
    // }

    logger().info("fast_validity_check done");

    // Reset the vertex and edge indices for the next patch
    for (auto vertex : vertices) {
      vertex_to_idx[vertex] = -1;
    }
    for (auto edge : edges) {
      edge_to_idx[edge] = -1;
    }
  }
}
