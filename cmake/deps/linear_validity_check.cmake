## ================== [linear_validity_check] ============
set(LINEAR_VALIDITY_CHECK_SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/linear_validity_check")

# Find system CGAL (newer version)
find_package(CGAL REQUIRED COMPONENTS Core)

add_library(linear_validity_check STATIC)

# Automatically include all cpp files in the directory
file(GLOB LINEAR_VALIDITY_CHECK_SOURCES
  "${LINEAR_VALIDITY_CHECK_SOURCE_DIR}/*.cpp"
)

target_sources(linear_validity_check PRIVATE
  ${LINEAR_VALIDITY_CHECK_SOURCES}
)
target_link_libraries(linear_validity_check PUBLIC
  CGAL::CGAL
  CGAL::CGAL_Core
  Eigen3::Eigen
  igl::igl
)
target_include_directories(linear_validity_check PUBLIC
  ${LINEAR_VALIDITY_CHECK_SOURCE_DIR}
)
add_library(linear_validity_check::linear_validity_check ALIAS linear_validity_check)
set_target_properties(linear_validity_check PROPERTIES FOLDER third_party)
