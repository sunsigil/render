#ifndef GEO_H
#define GEO_H

#include "linalg.h"

mat_t barycentric(mat_t tri, mat_t p);
mat_t inverse_barycentric(mat_t tri, mat_t p);
mat_t centroid(mat_t tri);
mat_t tri_AABB(mat_t tri);
int AABB_intersect(mat_t aabb_1, mat_t aabb_2);
double point_line_dist(mat_t p, mat_t a, mat_t b);
mat_t homogenize(mat_t tri);
int tri_contains(mat_t tri, mat_t p);
int AABB_contains(mat_t aabb, mat_t p);
mat_t mesh_AABB(double* v, int n);

#endif
