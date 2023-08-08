#include "draw.h"
#include "geo.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <float.h>

double sat(double n)
{
	if(n < 0.0)
	{ return 0.0; }
	else if(n > 1.0)
	{ return 1.0; }
	return n;
}

mat_t clr2mat(colour_t c)
{
	double r = c.r / 255.0;
	double g = c.g / 255.0;
	double b = c.b / 255.0;

	return mat_row(3, sat(r), sat(g), sat(b));
}

colour_t mat2clr(mat_t m)
{
	unsigned char r = sat(m.data[0]) * 255;
	unsigned char g = sat(m.data[1]) * 255;
	unsigned char b = sat(m.data[2]) * 255;

	colour_t clr = {r, g, b, 255};
	return clr;
}

mat_t sample(BMP_t* tex, mat_t uv)
{
	int x = uv.data[0] * (tex->width-1);
	int y = uv.data[1] * (tex->height-1);
	return clr2mat(BMP_get_pixel(tex, x, y));
}

mat_t cube_sample(BMP_t* tex, mat_t xyz)
{
	double x = xyz.data[0];
	double y = -xyz.data[1];
	double z = xyz.data[2];
	double xa = fabs(x);
	double ya = fabs(y);
	double za = fabs(z);

	double xu = 0.5 - (z/x) * 0.5;
	double xv = 0.5 - (y/x) * (x > 0.0 ? 0.5 : -0.5);
	double yu = 0.5 + (x/y) * (y > 0.0 ? 0.5 : -0.5);
	double yv = 0.5 - (z/y) * 0.5;
	double zu = 0.5 + (x/z) * 0.5;
	double zv = 0.5 - (y/z) * (z > 0.0 ? 0.5 : -0.5);

	mat_t xuv = mat_row(2, xu * 0.25, xv / 3.0);
	mat_t yuv = mat_row(2, yu * 0.25, (1.0-yv) / 3.0);
	mat_t zuv = mat_row(2, zu * 0.25, zv / 3.0);
	xuv = mat_add(xuv, mat_row(2, x > 0.0 ? 0.5 : 0.0, 1.0/3.0));
	yuv = mat_add(yuv, mat_row(2, 0.25, y > 0.0 ? 0.0 : 2.0/3.0));
	zuv = mat_add(zuv, mat_row(2, z > 0.0 ? 0.25 : 0.75, 1.0/3.0));

	double maxa = fmax(za, fmax(ya, xa));
	if(maxa == xa)
	{ return sample(tex, xuv); }
	if(maxa == ya)
	{ return sample(tex, yuv); }
	return sample(tex, zuv);
}

void point(BMP_t* bmp, colour_t clr, int x, int y, double r)
{
	mat_t cen = mat_row(2, (double) x, (double) y);

	for(int row = y-r; row < y+r; row++)
	{
		if(row < 0 || row >= bmp->height)
		{ continue; }

		for(int col = x-r; col < x+r; col++)
		{
			if(col < 0 || col >= bmp->width)
			{ continue; }

			mat_t pt = mat_row(2, (double) col, (double) row);
			mat_t ray = mat_sub(pt, cen);
			double dist = mat_norm(ray);

			if(dist <= r)
			{ BMP_set_pixel(bmp, col, row, clr); }
		}
	}
}

// implementation based on Dmitri Sokolov's
void line(BMP_t* bmp, colour_t clr, int x0, int y0, int x1, int y1)
{
	// if the line is steep, transpose its start and end points
	bool steep = abs(y1-y0) > abs(x1-x0);
	if(steep)
	{
		int temp = x0;
		x0 = y0;
		y0 = temp;

		temp = x1;
		x1 = y1;
		y1 = temp;
	}

	// if the line heads left, swap its start and end points
	bool leftward = x0 > x1;
	if(leftward)
	{
		int temp = x0;
		x0 = x1;
		x1 = temp;

		temp = y0;
		y0 = y1;
		y1 = temp;
	}
	
	int dx = x1 - x0;
	int dy = y1 - y0;

	// account for line heading up or down
	int y_step = (y1 > y0) ? 1 : -1;
	int y = y0;
	
	// approximate d_err as abs(dy) / (dx ~= 0.5)
	int d_err = abs(dy) * 2;
	int err = 0;

	// if line is steep, we swap x,y in the draw call to undo our earlier transposition
	// we employ a branch between two for loops to avoid branching within one loop
	if(steep)
	{
		for(int x = x0; x < x1; x++)
		{
			if(y >= 0 && y < bmp->width && x >= 0 && x < bmp->height)
			{ BMP_set_pixel(bmp, y, x, clr); }

			err += d_err;
			if(err > dx)
			{
				y += y_step;
				err -= dx*2;
			}
		}
	}
	else
	{
		for(int x = x0; x < x1; x++)
		{
			if(x >= 0 && x < bmp->width && y >= 0 && y < bmp->height)
			{ BMP_set_pixel(bmp, x, y, clr); }

			err += d_err;
			if(err > dx)
			{
				y += y_step;
				err -= dx*2;
			}
		}
	}
}

void tri(BMP_t* bmp, colour_t clr, mat_t tri)
{
	mat_t a = mat_get_row(tri, 0);
	mat_t b = mat_get_row(tri, 1);
	mat_t c = mat_get_row(tri, 2);

	line(bmp, clr, a.data[0], a.data[1], b.data[0], b.data[1]);
	line(bmp, clr, b.data[0], b.data[1], c.data[0], c.data[1]);
	line(bmp, clr, c.data[0], c.data[1], a.data[0], a.data[1]);
}

void AABB(BMP_t* bmp, colour_t clr, mat_t aabb, mat_t G)
{
	mat_t center = mat_get_row(aabb, 0);
	mat_t radii = mat_get_row(aabb, 1);

	double mx = center.data[0] - radii.data[0];
	double Mx = center.data[0] + radii.data[0];
	double my = center.data[1] - radii.data[1];
	double My = center.data[1] + radii.data[1];
	double mz = center.data[2] - radii.data[2];
	double Mz = center.data[2] + radii.data[2];

	mat_t mmm = mat_row(4, mx, my, mz, 1.0);
	mat_t mmM = mat_row(4, mx, my, Mz, 1.0);
	mat_t mMm = mat_row(4, mx, My, mz, 1.0);
	mat_t mMM = mat_row(4, mx, My, Mz, 1.0);
	mat_t Mmm = mat_row(4, Mx, my, mz, 1.0);
	mat_t MmM = mat_row(4, Mx, my, Mz, 1.0);
	mat_t MMm = mat_row(4, Mx, My, mz, 1.0);
	mat_t MMM = mat_row(4, Mx, My, Mz, 1.0);
	mmm = mat_mul(mmm, mat_trans(G));
	mmM = mat_mul(mmM, mat_trans(G));
	mMm = mat_mul(mMm, mat_trans(G));
	mMM = mat_mul(mMM, mat_trans(G));
	Mmm = mat_mul(Mmm, mat_trans(G));
	MmM = mat_mul(MmM, mat_trans(G));
	MMm = mat_mul(MMm, mat_trans(G));
	MMM = mat_mul(MMM, mat_trans(G));
	mmm = homogenize(mmm);
	mmM = homogenize(mmM);
	mMm = homogenize(mMm);
	mMM = homogenize(mMM);
	Mmm = homogenize(Mmm);
	MmM = homogenize(MmM);
	MMm = homogenize(MMm);
	MMM = homogenize(MMM);

	line(bmp, clr, mmm.data[0], mmm.data[1], Mmm.data[0], Mmm.data[1]);
	line(bmp, clr, mmm.data[0], mmm.data[1], mMm.data[0], mMm.data[1]);
	line(bmp, clr, mmm.data[0], mmm.data[1], mmM.data[0], mmM.data[1]);
	line(bmp, clr, MMM.data[0], MMM.data[1], MMm.data[0], MMm.data[1]);
	line(bmp, clr, MMM.data[0], MMM.data[1], MmM.data[0], MmM.data[1]);
	line(bmp, clr, MMM.data[0], MMM.data[1], mMM.data[0], mMM.data[1]);
	line(bmp, clr, mMm.data[0], mMm.data[1], MMm.data[0], MMm.data[1]);
	line(bmp, clr, MmM.data[0], MmM.data[1], mmM.data[0], mmM.data[1]);
	line(bmp, clr, mMm.data[0], mMm.data[1], mMM.data[0], mMM.data[1]);
	line(bmp, clr, mMM.data[0], mMM.data[1], mmM.data[0], mmM.data[1]);
	line(bmp, clr, MmM.data[0], MmM.data[1], Mmm.data[0], Mmm.data[1]);
	line(bmp, clr, Mmm.data[0], Mmm.data[1], MMm.data[0], MMm.data[1]);
}

void gizmo(BMP_t* bmp, mat_t G)
{
	mat_t basis = mat_id(4);
	basis = mat_set_col(basis, 3, mat_col(4, 1.0, 1.0, 1.0, 1.0));
	basis = mat_mul(basis, mat_trans(G));

	mat_t u = mat_get_row(basis, 0);
	mat_t v = mat_get_row(basis, 1);
	mat_t w = mat_get_row(basis, 2);
	mat_t o = mat_get_row(basis, 3);
	u = homogenize(u);
	v = homogenize(v);
	w = homogenize(w);
	o = homogenize(o);

	colour_t r = {255, 0, 0, 255};
	colour_t g = {0, 255, 0, 255};
	colour_t b = {0, 0, 255, 255};

	line(bmp, r, o.data[0], o.data[1], u.data[0], u.data[1]);
	line(bmp, g, o.data[0], o.data[1], v.data[0], v.data[1]);
	line(bmp, b, o.data[0], o.data[1], w.data[0], w.data[1]);
}

void sky(screen_t screen, uniforms_t uni)
{
	mat_t Tinv = mat4_inv(uni.T);

	double device_z = -1.0;

	for(int y = 0; y < screen.height; y++)
	{
		double half_height = (double) screen.height * 0.5;
		double device_y = ((double) y - half_height) / half_height;

		for(int x = 0; x < screen.width; x++)
		{
			double half_width = (double) screen.width * 0.5;
			double device_x = ((double) x - half_width) / half_width;

			mat_t ndc = mat_row(4, device_x, device_y, device_z, 1.0);
			double clip_w = 1.0 / mat_mul(ndc, mat_trans(Tinv)).data[3];
			mat_t world = mat_scale(mat_mul(ndc, mat_trans(Tinv)), clip_w);

			mat_t ray = mat_sub(world, uni.eye);
			mat_t direct = mat_normalized(ray);
			colour_t clr = mat2clr(cube_sample(uni.sky_tex, direct));

			BMP_set_pixel(screen.frame_tex, x, y, clr);
		}
	}
}

vert_out_t vert(uniforms_t uni, vert_in_t in)
{
	vert_out_t out;

	mat_t pos = mat_mul(in.pos, mat_trans(uni.T));
	out.pos_camera = pos;
	pos = homogenize(pos);
	out.pos_clip = pos;
	pos = mat_mul(pos, mat_trans(uni.W));
	out.pos_screen = pos;
	out.uv = in.uv;
	out.norm = mat_mul(in.norm, mat_trans(uni.N));

	return out;
}

void frag(screen_t screen, uniforms_t uni, vert_out_t in)
{
	mat_t a = mat_get_row(in.pos_screen, 0);
	mat_t b = mat_get_row(in.pos_screen, 1);
	mat_t c = mat_get_row(in.pos_screen, 2);

	mat_t xs = mat_get_col(in.pos_screen, 0);
	mat_t ys = mat_get_col(in.pos_screen, 1);
	mat_t zs = mat_get_col(in.pos_screen, 2);

	double max_x = mat_max(xs);
	double min_x = mat_min(xs);
	double max_y = mat_max(ys);
	double min_y = mat_min(ys);
	double max_z = mat_max(zs);
	double min_z = mat_min(zs);

	for(int y = min_y; y < max_y; y++)
	{
		if(y < 0 || y >= screen.height)
		{ continue; }

		for(int x = min_x; x < max_x; x++)
		{
			if(x < 0 || x >= screen.width)
			{ continue; }

			mat_t pixel = mat_row(4, (double) x, (double) y, 0.0, 1.0);
			mat_t tri = mat_set_col(in.pos_screen, 2, mat_col(3, 0.0, 0.0, 0.0));
			if(!tri_contains(in.pos_screen, pixel))
			{ continue; }

			mat_t bary_screen = barycentric(in.pos_screen, pixel);
			if
			(
			 	mat_min(bary_screen) > 0.0 &&
				mat_max(bary_screen) < 1.0 &&
				mat_validate(bary_screen)
			)
			{
				mat_t zs_cam = mat_get_col(in.pos_camera, 2); 
				double zpc_inv =
					bary_screen.data[0] / zs_cam.data[0] +
					bary_screen.data[1] / zs_cam.data[1] +
					bary_screen.data[2] / zs_cam.data[2];
				double zpc = 1.0 / zpc_inv;
				double vpc = bary_screen.data[1] * zpc / zs_cam.data[1];
				double wpc = bary_screen.data[2] * zpc / zs_cam.data[2];
				double upc = 1.0 - vpc - wpc;
				mat_t bary = mat_row(3, upc, vpc, wpc);
				
				double z = mat_lincomb(zs, bary_screen).data[0];
				double buf_z = screen.depth_buffer[y * screen.width + x];

				if(z > buf_z)
				{
					mat_t pos = inverse_barycentric(in.pos_camera, bary_screen);
					mat_t uv = inverse_barycentric(in.uv, bary);
					mat_t norm = inverse_barycentric(in.norm, bary);
					
					norm = mat_normalized(norm);
					mat_t light = mat_normalized(mat_scale(uni.beam, -1.0));
					mat_t view = mat_normalized(mat_sub(uni.eye, pos));
					
					double nl = sat(mat_dot(norm, light));
					mat_t direct = mat_row(3, 1.0, 1.0, 1.0);
					mat_t ambient = mat_scale(direct, 0.15);
					
					double albedo = 1.0;
					mat_t diffuse = sample(uni.mesh_tex, uv);
			
					mat_t half = mat_normalized(mat_add(light, view));
					double nh = sat(mat_dot(norm, half));
					double highlight = nl > 0.0 ? pow(nh, 20.0) : 0.0;
					mat_t specular = direct;

					diffuse = mat_scale(diffuse, albedo / M_PI * nl);
					specular = mat_scale(specular, highlight);

					mat_t surface = mat_compwise(mat_add(diffuse, specular), direct);
					surface = mat_add(surface, ambient);

					BMP_set_pixel(screen.frame_tex, x, y, mat2clr(surface));
					BMP_set_pixel(screen.frame_tex, x, y, mat2clr(sample(uni.mesh_tex, uv)));
					screen.depth_buffer[y * screen.width + x] = z;
				}
			}
		}
	}
}

void depth(screen_t screen)
{
	double min = DBL_MAX;
	double max = -DBL_MAX;
	for(int i = 0; i < screen.width*screen.height; i++)
	{
		min = fmin(min, screen.depth_buffer[i]);
		max = fmax(max, screen.depth_buffer[i]);
	}

	for(int y = 0; y < screen.height; y++)
	{
		for(int x = 0; x < screen.width; x++)
		{
			double d = screen.depth_buffer[y * screen.width + x];
			d = (d - min) / (max - min);

			mat_t mat = mat_row(4, d, d, d, 1.0);
			BMP_set_pixel(screen.depth_tex, x, y, mat2clr(mat));
		}
	}
}

