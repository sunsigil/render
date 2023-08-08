#ifndef DRAW_H
#define DRAW_H

#include "BMP.h"
#include "linalg.h"

typedef struct screen
{
	BMP_t* frame_tex;
	double* depth_buffer;
	BMP_t* depth_tex;

	int width;
	int height;
} screen_t;

typedef struct uniforms
{
	BMP_t* sky_tex;
	BMP_t* mesh_tex;

	mat_t eye;
	mat_t look;
	mat_t beam;

	mat_t T;
	mat_t W;
	mat_t N;
} uniforms_t;

typedef struct vert_in
{
	mat_t pos;
	mat_t uv;
	mat_t norm;
} vert_in_t;

typedef struct vert_out
{
	mat_t pos_camera;
	mat_t pos_clip;
	mat_t pos_screen;
	mat_t uv;
	mat_t norm;
} vert_out_t;

void point(BMP_t* bmp, colour_t clr, int x, int y, double r);
void line(BMP_t* bmp, colour_t clr, int x0, int y0, int x1, int y1);
void tri(BMP_t* bmp, colour_t clr, mat_t tri);
void AABB(BMP_t* bmp, colour_t clr, mat_t aabb, mat_t G);
void gizmo(BMP_t* bmp, mat_t G);

void sky(screen_t screen, uniforms_t uni);
vert_out_t vert(uniforms_t uni, vert_in_t in);
void frag(screen_t screen, uniforms_t uni, vert_out_t in);
void depth(screen_t screen);

#endif
