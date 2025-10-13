// 1. triangulate each face w/ winding number > 0 in the arrangement
// 2. connect the vertices with the connected components
// 3. iteratively remove the internal vertices via edge flips (decimate_internal_vertices.h)
// 4. divide it into simple polygons using Shor van Wyck's `simplify_triangulation` function (simplify_polygon.h)
