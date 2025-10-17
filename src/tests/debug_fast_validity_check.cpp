#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <filesystem>

#include "linear_validity_check/fast_validity_check.h"

// Function to read Eigen matrix from file
template<typename MatrixType>
MatrixType readMatrixFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }
    
    std::vector<std::vector<double>> data;
    std::string line;
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        std::vector<double> row;
        std::istringstream iss(line);
        double value;
        while (iss >> value) {
            row.push_back(value);
        }
        if (!row.empty()) {
            data.push_back(row);
        }
    }
    
    if (data.empty()) {
        throw std::runtime_error("No data found in file: " + filename);
    }
    
    const int rows = data.size();
    const int cols = data[0].size();
    
    // Check all rows have same number of columns
    for (const auto& row : data) {
        if (static_cast<int>(row.size()) != cols) {
            throw std::runtime_error("Inconsistent column count in file: " + filename);
        }
    }
    
    MatrixType matrix(rows, cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            matrix(i, j) = data[i][j];
        }
    }
    
    return matrix;
}

// Function to read Eigen vector from file
template<typename VectorType>
VectorType readVectorFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }
    
    std::vector<double> data;
    std::string line;
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        std::istringstream iss(line);
        double value;
        while (iss >> value) {
            data.push_back(value);
        }
    }
    
    if (data.empty()) {
        throw std::runtime_error("No data found in file: " + filename);
    }
    
    VectorType vector(static_cast<int>(data.size()));
    for (int i = 0; i < static_cast<int>(data.size()); i++) {
        vector(i) = data[i];
    }
    
    return vector;
}

// Function to extract ID from debug directory path
std::string extractIdFromPath(const std::string& debug_dir) {
    std::filesystem::path path(debug_dir);
    std::string dir_name = path.filename().string();
    
    // Remove trailing slash if present
    if (!dir_name.empty() && dir_name.back() == '/') {
        dir_name.pop_back();
    }
    
    return dir_name;
}

// Function to export mesh as OBJ file
void exportMeshAsObj(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not create file: " + filename);
    }
    
    // Write vertices
    for (int i = 0; i < V.rows(); i++) {
        file << "v " << V(i, 0) << " " << V(i, 1) << " " << 0.000 << std::endl;
    }
    
    // Write faces (OBJ uses 1-based indexing)
    for (int i = 0; i < F.rows(); i++) {
        file << "f";
        for (int j = 0; j < F.cols(); j++) {
            file << " " << (F(i, j) + 1); // Convert to 1-based indexing
        }
        file << std::endl;
    }
    
    file.close();
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <debug_inputs_directory>" << std::endl;
        std::cerr << "Example: " << argv[0] << " debug_inputs" << std::endl;
        return 1;
    }
    
    std::string debug_dir = argv[1];
    
    try {
        std::cout << "Loading debug inputs from directory: " << debug_dir << std::endl;
        
        // Read all the input matrices and vectors
        auto vertices_3d = readMatrixFromFile<Eigen::MatrixXd>(debug_dir + "/vertices_3d.txt");
        auto vertices_2d = readMatrixFromFile<Eigen::MatrixXd>(debug_dir + "/vertices_2d.txt");
        auto edges_connectivity = readMatrixFromFile<Eigen::MatrixXi>(debug_dir + "/edges_connectivity.txt");
        auto edges_sign = readVectorFromFile<Eigen::VectorXi>(debug_dir + "/edges_sign.txt");
        auto edges_is_cut = readVectorFromFile<Eigen::VectorXi>(debug_dir + "/edges_is_cut.txt");
        auto vertex_is_cusp = readVectorFromFile<Eigen::VectorXi>(debug_dir + "/vertex_is_cusp.txt");
        auto vertex_is_singularity = readVectorFromFile<Eigen::VectorXi>(debug_dir + "/vertex_is_singularity.txt");
        auto camera_pos = readVectorFromFile<Eigen::Vector3d>(debug_dir + "/camera_pos.txt");
        
        std::cout << "Successfully loaded all input data:" << std::endl;
        std::cout << "  vertices_3d: " << vertices_3d.rows() << "x" << vertices_3d.cols() << std::endl;
        std::cout << "  vertices_2d: " << vertices_2d.rows() << "x" << vertices_2d.cols() << std::endl;
        std::cout << "  edges_connectivity: " << edges_connectivity.rows() << "x" << edges_connectivity.cols() << std::endl;
        std::cout << "  edges_sign: " << edges_sign.size() << std::endl;
        std::cout << "  edges_is_cut: " << edges_is_cut.size() << std::endl;
        std::cout << "  vertex_is_cusp: " << vertex_is_cusp.size() << std::endl;
        std::cout << "  vertex_is_singularity: " << vertex_is_singularity.size() << std::endl;
        std::cout << "  camera_pos: " << camera_pos.transpose() << std::endl;
        
        std::cout << "\nCalling fast_validity_check..." << std::endl;
        
        // Call the fast_validity_check function
        auto [is_wso_succeeded, V_out, F_out, V_JI] = utils::fast_validity_check(
            vertices_3d, vertices_2d, edges_connectivity, edges_sign, edges_is_cut,
            vertex_is_cusp, vertex_is_singularity, camera_pos, true);

        std::cout << "fast_validity_check completed!" << std::endl;
        std::cout << "Result: " << (is_wso_succeeded ? "SUCCESS" : "FAILED")
                  << std::endl;

        // export the mesh
        std::string id = extractIdFromPath(debug_dir);
        std::string output_dir = "./debug_output/" + id;
        
        // Create output directory if it doesn't exist
        std::filesystem::create_directories(output_dir);
        
        if (is_wso_succeeded) {
            // Export mesh as OBJ file
            std::string mesh_filename = output_dir + "/mesh.obj";
            exportMeshAsObj(V_out, F_out, mesh_filename);
            std::cout << "Mesh exported to: " << mesh_filename << std::endl;
        }
        
        // Move SVG file if it exists
        std::string svg_source = "./.tmp_debug/validity_check.svg";
        std::string svg_dest = output_dir + "/validity_check.svg";
        
        if (std::filesystem::exists(svg_source)) {
            std::filesystem::rename(svg_source, svg_dest);
            std::cout << "SVG file moved to: " << svg_dest << std::endl;
        } else {
            std::cout << "Warning: SVG file not found at " << svg_source << std::endl;
        }
        
        
        if (is_wso_succeeded) {
            std::cout << "  V_out: " << V_out.rows() << "x" << V_out.cols() << std::endl;
            std::cout << "  F_out: " << F_out.rows() << "x" << F_out.cols() << std::endl;
            std::cout << "  V_JI: " << V_JI.size() << std::endl;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
