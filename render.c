#include "linalg.h"
#include "draw.h"
#include "geo.h"
#include "obj.h"
#include "vbo.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include <string.h>

#define DEG2RAD 0.01745329252

#define WIDTH 1920
#define HEIGHT 1080
#define MESH_PATH "mesh/currybox.obj"
#define TEX_PATH "img/tex/gcbox.bmp"
#define SKY_PATH "img/tex/checker.bmp"

double frand(float min, float max)
{
	double core = (double) rand() / (double) RAND_MAX;
	double range = max - min;
	return min + core * range;
}

void startwatch(clock_t* watch)
{ *watch = clock(); }
double stopwatch(clock_t* watch)
{ return (double) (clock() - *watch) / CLOCKS_PER_SEC; }

int main()
{
	colour_t red = { 255, 0, 0, 255 };
	colour_t green = { 0, 255, 0, 255 };
	colour_t blue = { 0, 0, 255, 255 };
	colour_t white = { 255, 255, 255, 255 };

	mat_t right = mat_get_row(mat_id(4), 0);
	mat_t up = mat_get_row(mat_id(4), 1);
	mat_t forward = mat_get_row(mat_id(4), 2);

	mat_t origin = mat_get_row(mat_id(4), 3);
	mat_t light = mat_row(4, 0.0, -1.0, -1.0, 0.0);
	mat_t eye = mat_row(4, 10.0, 0.0, -10.0, 1.0);

	double aspect = (double) WIDTH / (double) HEIGHT;
	double fovy = 60 * DEG2RAD;
	double near = -0.01;
	double far = -1000.0;
	
	mat_t C = lookat(eye, origin, up);
	mat_t P = perspective(aspect, fovy, near, far);
	mat_t T = mat_mul(P, C);
	mat_t N = mat_blit(C, 0, 3, mat_init(3, 1));
	N = mat_trans(mat4_inv(N));
	mat_t W = window(WIDTH, HEIGHT);
	mat_t G = mat_mul(W, T);

	time_t t = time(NULL);
	srand((unsigned int) t);
	clock_t watch;

	startwatch(&watch);	
	OBJ_t obj = OBJ_read(MESH_PATH);
	VBO_t vbo = VBO_init(obj);
	OBJ_dispose(&obj);
	
	BMP_t frame_tex = BMP_create(WIDTH, HEIGHT, 3);
	double* depth_buffer = calloc(WIDTH * HEIGHT, sizeof(double));
	BMP_t depth_tex = BMP_create(WIDTH, HEIGHT, 3);

	BMP_t mesh_tex = BMP_read(TEX_PATH); 
	BMP_t sky_tex = BMP_read(SKY_PATH);
	printf("LOAD %fs\n", stopwatch(&watch));

	screen_t screen;
	screen.frame_tex = &frame_tex;
	screen.depth_buffer = depth_buffer;
	screen.depth_tex = &depth_tex;
	screen.width = frame_tex.width;
	screen.height = frame_tex.height;

	uniforms_t uni;
	uni.sky_tex = &sky_tex;
	uni.mesh_tex = &mesh_tex;
	uni.eye = eye;
	uni.look = mat_normalized(mat_sub(origin, eye));
	uni.beam = light;
	uni.T = T;
	uni.W = W;
	uni.N = N;

	double min_time = DBL_MAX;
	double max_time = -DBL_MAX;
	double draw_time = 0.0;

	startwatch(&watch);

	sky(screen, uni);
	for(int i = 0; i < vbo.face_count; i++)
	{
		int fidx = i*3;
		
		vert_in_t in;
		in.pos = mat_init(3, 4);
		in.uv = mat_init(3, 2);
		in.norm = mat_init(3, 4);

		for(int pt = 0; pt < 3; pt++)
		{
			int vidx = vbo.indices[fidx+pt]*8;

			double x = vbo.vertices[vidx+0];
			double y = vbo.vertices[vidx+1];
			double z = vbo.vertices[vidx+2];

			double u = vbo.vertices[vidx+3];
			double v = vbo.vertices[vidx+4];

			double a = vbo.vertices[vidx+5];
			double b = vbo.vertices[vidx+6];
			double c = vbo.vertices[vidx+7];
			
			in.pos = mat_set_row(in.pos, pt, mat_row(4, x, y, z, 1.0));
			in.uv = mat_set_row(in.uv, pt, mat_row(2, u, v));
			in.norm = mat_set_row(in.norm, pt, mat_row(4, a, b, c, 0.0));
		}

		vert_out_t out = vert(uni, in);
		
		bool clip = false;
		for(int i = 0; i < 3; i++)
		{
			mat_t pt_clip = mat_get_row(out.pos_clip, i);
			if(pt_clip.data[0] < -1.0 || pt_clip.data[0] > 1.0)
			{ clip = true; break; }
			if(pt_clip.data[1] < -1.0 || pt_clip.data[1] > 1.0)
			{ clip = true; break; }
			if(pt_clip.data[2] < -1.0 || pt_clip.data[2] > 1.0)
			{ clip = true; break; }
		}
		if(clip)
		{ continue; }

		mat_t norm = mat_get_row(out.norm, 0);
		mat_t cen = centroid(out.pos_camera);
		mat_t view = mat_scale(cen, -1.0);
		double align = mat_dot(view, mat_get_row(norm, 0));
		if(align < 0.0)
		{ continue; }

		frag(screen, uni, out); 
		// tri(&frame_tex, white, out.pos_screen);

		/*cen = mat_mul(cen, mat_trans(W));
		cen = homogenize(cen);
		norm = mat_scale(norm, HEIGHT / 100.0);
		norm = mat_add(cen, norm);
		line(&frame_tex, red,  cen.data[0], cen.data[1], norm.data[0], norm.data[1]);*/
	}
	// gizmo(&frame_tex, G);
	// AABB(&frame_tex, white, mesh_AABB(obj.vs, obj.v_count), G);
	depth(screen);

	printf("DRAW %fs\n", stopwatch(&watch));

	BMP_write(&frame_tex, "img/out/frame_tex.bmp");
	BMP_write(&depth_tex, "img/out/depth_tex.bmp"); 
	
	free(depth_buffer);
	BMP_dispose(&depth_tex);
	BMP_dispose(&frame_tex);
	BMP_dispose(&mesh_tex);
	BMP_dispose(&sky_tex);

	return 0;
}
