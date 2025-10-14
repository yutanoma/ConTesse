#pragma once

#include "mesh.h"
#include <unordered_set>

void linear_wso_check(
  Mesh const &mesh, Camera const &camera,
  std::map<size_t, Mesh> &patch_triangulations, bool do_refine = false,
  std::unordered_set<int> selected_patches = std::unordered_set<int>(),
  bool ff_only = false);
