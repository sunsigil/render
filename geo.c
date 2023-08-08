#include "geo.h"
#include <stdlib.h>
#include <math.h>
#include <float.h>

mat_t barycentric(mat_t tri, mat_t p)
{
	// Real-Time Collision Detection, Christer Ericson
	
	double* a = tri.data;
	double* b = a+4;
	double* c = b+4;
	double* pp = p.data;

	double a0 = a[0]; double a1 = a[1]; double a2 = a[2];
	double b0 = b[0]; double b1 = b[1]; double b2 = b[2];
	double c0 = c[0]; double c1 = c[1]; double c2 = c[2];
	double p0 = pp[0]; double p1 = pp[1]; double p2 = pp[2];

	double v00 = b0-a0; double v01 = b1-a1; double v02 = b2-a2;
	double v10 = c0-a0; double v11 = c1-a1; double v12 = c2-a2;
	double v20 = p0-a0; double v21 = p1-a1; double v22 = p2-a2;

	double d00 = v00*v00 + v01*v01 + v02*v02; 
	double d01 = v00*v10 + v01*v11 + v02*v12; 
	double d11 = v10*v10 + v11*v11 + v12*v12; 
	double d20 = v20*v00 + v21*v01 + v22*v02; 
	double d21 = v20*v10 + v21*v11 + v22*v12; 
	double denom = (d00 * d11) - (d01 * d01);

	double v = ((d11 * d20) - (d01 * d21)) / denom;
	double w = ((d00 * d21) - (d01 * d20)) / denom;
	double u = 1.0 - v - w;

	return mat_row(3, u, v, w);
}

mat_t inverse_barycentric(mat_t tri, mat_t p)
{
	mat_t a = mat_get_row(tri, 0);
	mat_t b = mat_get_row(tri, 1);
	mat_t c = mat_get_row(tri, 2);

	mat_t ua = mat_scale(a, p.data[0]);
	mat_t vb = mat_scale(b, p.data[1]);
	mat_t wc = mat_scale(c, p.data[2]);

	return mat_add(mat_add(ua, vb), wc);
}

mat_t centroid(mat_t tri)
{
	double third = 1.0 / 3.0;
	mat_t thirds = mat_row(3, third, third, third);
	return inverse_barycentric(tri, thirds);
}

mat_t tri_AABB(mat_t tri)
{
	mat_t a = mat_get_row(tri, 0);
	mat_t b = mat_get_row(tri, 1);
	mat_t c = mat_get_row(tri, 2);
	
	double min_x = fmin(fmin(a.data[0], b.data[0]), c.data[0]);
	double min_y = fmin(fmin(a.data[1], b.data[1]), c.data[1]);
	double min_z = fmin(fmin(a.data[2], b.data[2]), c.data[2]);
	double max_x = fmax(fmax(a.data[0], b.data[0]), c.data[0]);
	double max_y = fmax(fmax(a.data[1], b.data[1]), c.data[1]);
	double max_z = fmax(fmax(a.data[2], b.data[2]), c.data[2]);

	double radius_x = (max_x - min_x) * 0.5;
	double radius_y = (max_y - min_y) * 0.5;
	double radius_z = (max_z - min_z) * 0.5;
	double center_x = min_x + radius_x;
	double center_y = min_y + radius_y;
	double center_z = min_z + radius_z;
	mat_t center = mat_row(3, center_x, center_y, center_z);
	mat_t radii = mat_row(3, radius_x, radius_y, radius_z);

	mat_t aabb = mat_init(2, 3);
	aabb = mat_set_row(aabb, 0, center);
	aabb = mat_set_row(aabb, 1, radii);
	
	return aabb;
}

int AABB_intersect(mat_t aabb_1, mat_t aabb_2)
{
	mat_t c_1 = mat_get_row(aabb_1, 0);
	mat_t c_2 = mat_get_row(aabb_2, 0);
	mat_t r_1 = mat_get_row(aabb_1, 1);
	mat_t r_2 = mat_get_row(aabb_2, 1);

	if(abs(c_1.data[0] - c_2.data[0]) > (r_1.data[0] + r_2.data[0]))
	{ return 0; }
	if(abs(c_1.data[1] - c_2.data[1]) > (r_1.data[1] + r_2.data[1]))
	{ return 0; }
	if(abs(c_1.data[2] - c_2.data[2]) > (r_1.data[2] + r_2.data[2]))
	{ return 0; }
	
	return 1;
}

double point_line_dist(mat_t p, mat_t a, mat_t b)
{
	mat_t ap = mat_sub(p, a);
	mat_t ab = mat_sub(b, a);

	double apap = mat_dot(ap, ap);
	double apab = mat_dot(ap, ab);
	double abab = mat_dot(ab, ab);

	return  sqrt(apap - apab * apab / abab);
}

mat_t homogenize(mat_t tri)
{
	for(int i = 0; i < 3; i++)
	{
		double* row = tri.data + i*4;
		double w_inv = 1.0 / row[3];
		row[0] *= w_inv;
		row[1] *= w_inv;
		row[2] *= w_inv;
		row[3] *= w_inv;
	}

	return tri;
}

int tri_contains(mat_t tri, mat_t p)
{
	mat_t a = mat_get_row(tri, 0);
	mat_t b = mat_get_row(tri, 1);	
	mat_t c = mat_get_row(tri, 2);	

	mat_t ab = mat_sub(b, a);
	mat_t ap = mat_sub(p, a);
	mat_t abap = vec3_cross(ab, ap);
	if(abap.data[2] < 0.0)
	{ return 0; }

	mat_t bc = mat_sub(c, b);
	mat_t bp = mat_sub(p, b);
	mat_t bcbp = vec3_cross(bc, bp);
	if(bcbp.data[2] < 0.0)
	{ return 0; }

	mat_t ca = mat_sub(a, c);
	mat_t cp = mat_sub(p, c);
	mat_t cacp = vec3_cross(ca, cp);
	if(cacp.data[2] < 0.0)
	{ return 0; }

	return 1;
}

int AABB_contains(mat_t aabb, mat_t p)
{
	mat_t center = mat_get_row(aabb, 0);
	mat_t radii = mat_get_row(aabb, 1);

	if(abs(p.data[0] - center.data[0]) > radii.data[0])
	{ return 0; }
	if(abs(p.data[1] - center.data[1]) > radii.data[1])
	{ return 0; }

	return 1;
}

mat_t mesh_AABB(double* v, int n)
{
	double min_x = DBL_MAX;
	double max_x = -DBL_MAX;
	double min_y = DBL_MAX;
	double max_y = -DBL_MAX;
	double min_z = DBL_MAX;
	double max_z = -DBL_MAX;
	
	for(int i = 0; i < n*3; i+=3)
	{
		min_x = fmin(min_x, v[i+0]);
		max_x = fmax(max_x, v[i+0]);
		min_y = fmin(min_y, v[i+1]);
		max_y = fmax(max_y, v[i+1]);
		min_z = fmin(min_z, v[i+2]);
		max_z = fmax(max_z, v[i+2]);
	}

	double cen_x = (min_x + max_x) * 0.5;
	double cen_y = (min_y + max_y) * 0.5;
	double cen_z = (min_z + max_z) * 0.5;
	double rad_x = max_x - cen_x;
	double rad_y = max_y - cen_y;
	double rad_z = max_z - cen_z;

	mat_t aabb = mat_init(2, 3);
	aabb = mat_set_row(aabb, 0, mat_row(3, cen_x, cen_y, cen_z));
	aabb = mat_set_row(aabb, 1, mat_row(3, rad_x, rad_y, rad_z));
	
	return aabb;
}
