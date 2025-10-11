if(TARGET igl::igl)
  return()
endif()

message(STATUS "contour-tesselation: Adding target igl::igl")

include(FetchContent)

## ================ IGL

FetchContent_Declare(
  libigl
  GIT_REPOSITORY https://github.com/libigl/libigl.git
  GIT_TAG c7324d3fa6fa94e55764b274ab15949ebbd2471f
)
FetchContent_GetProperties(libigl)
if(NOT libigl_POPULATED)
  FetchContent_Populate(libigl)
endif()

# Configure libigl options before including it
set(LIBIGL_WITH_CGAL ON CACHE BOOL "Use CGAL" FORCE)
set(LIBIGL_WITH_PREDICATES ON CACHE BOOL "Use predicates" FORCE)
set(LIBIGL_WITH_TRIANGLE ON CACHE BOOL "Use triangle" FORCE)
set(LIBIGL_WITH_OPENGL ON CACHE BOOL "Use OpenGL" FORCE)
set(LIBIGL_WITH_GLFW ON CACHE BOOL "Use GLFW" FORCE)
set(LIBIGL_WITH_EMBREE OFF CACHE BOOL "Use Embree" FORCE)
set(LIBIGL_WITH_IMGUI OFF CACHE BOOL "Use ImGui" FORCE)
set(LIBIGL_WITH_PNG OFF CACHE BOOL "Use PNG" FORCE)
set(LIBIGL_WITH_XML OFF CACHE BOOL "Use XML" FORCE)
set(LIBIGL_WITH_FWINDING_NUMBER OFF CACHE BOOL "Use Fast Winding Number" FORCE)

# Add libigl as a subdirectory to build it properly
add_subdirectory(${libigl_SOURCE_DIR} libigl_build EXCLUDE_FROM_ALL)

# Exclude problematic files for ARM architecture
if(APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "arm64")
    # Remove the problematic file from the build
    file(REMOVE ${libigl_SOURCE_DIR}/include/igl/fast_winding_number.cpp)
endif()

# Create an alias for the main igl target that includes all components
add_library(igl::igl ALIAS igl)

# The libigl subdirectory will create all the necessary targets:
# - igl::core
# - igl::predicates  
# - igl::copyleft::cgal
# - igl::restricted::triangle
# - igl::opengl
# - etc.

# kill warnings
# if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang" OR
#   "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang"
# )
#   target_compile_options(igl PRIVATE "-Wno-unused-function")
# endif()

## ================ igl::igl

add_library(igl::igl ALIAS igl)