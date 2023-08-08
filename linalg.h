#ifndef LINALG_H
#define LINALG_H

#define MAT_MAX_SIZE 16 
#define LINALG_EPS DBL_EPSILON
// #define LINALG_EPS 1e-15

typedef struct mat
{
	int m;
	int n;

	double data[MAT_MAX_SIZE];
} mat_t;

mat_t mat_init(int m, int n);
mat_t mat_col(int m, ...);
mat_t mat_row(int n, ...);
mat_t mat_id(int m);

int mat_validate(mat_t A);
mat_t mat_sanitize(mat_t A);
double mat_min(mat_t A);
double mat_max(mat_t A);
double mat_sum(mat_t A);
double mat_prod(mat_t A);
double mat_norm(mat_t A);
mat_t mat_normalized(mat_t A);
mat_t mat_trans(mat_t A);

mat_t mat_set(mat_t A, int i, int j, double x);
double mat_get(mat_t A, int i, int j);
mat_t mat_set_row(mat_t A, int i, mat_t R);
mat_t mat_set_col(mat_t A, int j, mat_t C);
mat_t mat_get_row(mat_t A, int i);
mat_t mat_get_col(mat_t A, int j);
mat_t mat_blit(mat_t A, int i, int j, mat_t B);

mat_t mat_reshape(mat_t A, int m, int n);
mat_t mat_resize(mat_t A, int m, int n);
mat_t mat_add_row(mat_t A, int i, mat_t R);
mat_t mat_add_col(mat_t A, int j, mat_t C);
mat_t mat_del_row(mat_t A, int i);
mat_t mat_del_col(mat_t A, int j);
mat_t mat_swap_rows(mat_t A, int i, int j);
mat_t mat_swap_cols(mat_t A, int i, int j);

mat_t mat_scale(mat_t A, double l);
mat_t mat_add(mat_t A, mat_t B);
mat_t mat_sub(mat_t A, mat_t B);
mat_t mat_mul(mat_t A, mat_t B);
double mat_dot(mat_t A, mat_t B);
mat_t mat_lerp(mat_t A, mat_t B, double t);
mat_t mat_compwise(mat_t A, mat_t B);
mat_t mat_lincomb(mat_t A, mat_t B);

double mat2_det(mat_t A);
mat_t vec3_cross(mat_t A, mat_t B);
mat_t mat4_inv(mat_t A);

mat_t scaling(double x, double y, double z);
mat_t rotation(double x, double y, double z);
mat_t translation(double x, double y, double z);
mat_t lookat(mat_t eye, mat_t target, mat_t up);
mat_t perspective(double aspect, double fov, double near, double far);
mat_t window(int width, int height);

void mat_print(mat_t A);
void mat_dump(mat_t A);
void num_print(double n);

#endif
