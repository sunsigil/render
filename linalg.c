#include "linalg.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <string.h>

mat_t mat_init(int m, int n)
{
	if(m*n > MAT_MAX_SIZE)
	{
		printf("[mat_init] product of matrix dimensions may not exceed %d\n", MAT_MAX_SIZE);
		exit(EXIT_FAILURE);
	}

	mat_t A;
	A.m = m;
	A.n = n;
	memset(A.data, 0, MAT_MAX_SIZE * sizeof(double));

	return A;
}

mat_t mat_row(int n, ...)
{
	mat_t R = mat_init(1, n);

	va_list args;
	va_start(args, n);
	for(int j = 0; j < n; j++)
	{
		R.data[j] = va_arg(args, double);
	}
	va_end(args);

	return R;
}

mat_t mat_col(int m, ...)
{
	mat_t C = mat_init(m, 1);
	
	va_list args;
	va_start(args, m);
	for(int i = 0; i < m; i++)
	{
		C.data[i] = va_arg(args, double);
	}
	va_end(args);

	return C;
}

mat_t mat_id(int m)
{
	mat_t A = mat_init(m, m);

	for(int i = 0; i < m*m; i+=m+1)
	{
		A.data[i] = 1;
	}

	return A;
}

int mat_validate(mat_t A)
{
	for(int i = 0; i < A.m * A.n; i++)
	{
		if(!isfinite(A.data[i]))
		{ return 0; }
	}
	return 1;
}

mat_t mat_sanitize(mat_t A)
{
	for(int i = 0; i < A.m * A.n; i++)
	{
		if(fabs(A.data[i]) < LINALG_EPS)
		{ A.data[i] = 0.0; }
	}
	return A;
}

double mat_min(mat_t A)
{
	double min = DBL_MAX;
	for(int i = 0; i < A.m * A.n; i++)
	{ min = min < A.data[i] ? min : A.data[i]; }
	return min;
}

double mat_max(mat_t A)
{
	double max = -DBL_MAX;
	for(int i = 0; i < A.m * A.n; i++)
	{ max = max > A.data[i] ? max : A.data[i]; }
	return max;
}

double mat_sum(mat_t A)
{
	double sum = 0.0;
	for(int i = 0; i < A.m * A.n; i++)
	{ sum += A.data[i]; }
	return sum;
}

double mat_prod(mat_t A)
{
	double prod = 1.0;
	for(int i = 0; i < A.m * A.n; i++)
	{ prod *= A.data[i]; }
	return prod;
}

mat_t mat_set(mat_t A, int i, int j, double x)
{
	if(i >= A.m || j >= A.n)
	{
		puts("[mat_set] out of bounds");
		exit(EXIT_FAILURE);
	}

	A.data[i * A.n + j] = x;

	return A;
}

double mat_get(mat_t A, int i, int j)
{
	if(i >= A.m || j >= A.n)
	{
		puts("[mat_set] index out of bounds");
		exit(EXIT_FAILURE);
	}

	return A.data[i * A.n + j];
}

mat_t mat_set_row(mat_t A, int i, mat_t R)
{
	if(i >= A.m)
	{
		puts("[mat_set_row] index out of bounds");
		exit(EXIT_FAILURE);
	}
	if(R.m != 1 || R.n != A.n)
	{
		puts("[mat_set_row] incorrect operand dimensions");
		exit(EXIT_FAILURE);
	}

	for(int j = 0; j < A.n; j++)
	{
		A.data[i * A.n + j] = R.data[j];
	}
	
	return A;
}

mat_t mat_set_col(mat_t A, int j, mat_t C)
{
	if(j >= A.n)
	{
		puts("[mat_set_col] index out of bounds");
		exit(EXIT_FAILURE);
	}
	if(C.m != A.m || C.n != 1)
	{
		puts("[mat_set_col] incorrect operand dimensions");
		exit(EXIT_FAILURE);
	}

	for(int i = 0; i < A.m; i++)
	{
		A.data[i * A.n + j] = C.data[i];
	}
	
	return A;
}

mat_t mat_get_row(mat_t A, int i)
{
	if(i >= A.m)
	{
		puts("[mat_get_row] index out of bounds");
		exit(EXIT_FAILURE);
	}

	mat_t R = mat_init(1, A.n);

	for(int j = 0; j < A.n; j++)
	{
		R.data[j] = A.data[i * A.n + j];
	}
	
	return R;
}

mat_t mat_get_col(mat_t A, int j)
{
	if(j < 0 || j >= A.n)
	{
		puts("[mat_get_col] index out of bounds");
		exit(EXIT_FAILURE);
	}

	mat_t C = mat_init(A.m, 1);

	for(int i = 0; i < A.m; i++)
	{
		C.data[i] = A.data[i * A.n + j];
	}

	return C;
}

mat_t mat_blit(mat_t A, int i, int j, mat_t B)
{
	if(i < 0 || j < 0 || i >= A.m || j >= A.n)
	{
		puts("[mat_blit] index out of bounds");
		exit(EXIT_FAILURE);
	}

	if(i + B.m > A.m || j + B.n > A.n)
	{
		puts("[mat_blit] blit oversteps bounds");
		exit(EXIT_FAILURE);
	}

	for(int row = i; row < i + B.m; row++)
	{
		for(int col = j; col < j + B.n; col++)
		{
			int b_row = row - i;
			int b_col = col - j;
			A.data[row * A.n + col] =
			B.data[b_row * B.n + b_col];
		}
	}
	
	return A;
}

mat_t mat_reshape(mat_t A, int m, int n)
{
	if(m*n != A.m * A.n)
	{
		printf("[mat_reshape] product of new matrix dimensions must equal product of old matrix dimensions %d\n", MAT_MAX_SIZE);
		exit(EXIT_FAILURE);
	}

	A.m = m;
	A.n = n;

	return A;
}

mat_t mat_resize(mat_t A, int m, int n)
{
	if(m < 0 || n < 0 || m*n > MAT_MAX_SIZE)
	{
		printf("[mat_resize] product of new matrix dimensions may not exceed %d\n", MAT_MAX_SIZE);
		exit(EXIT_FAILURE);
	}

	mat_t B = mat_init(m, n);
	for(int row = 0; row < fmin(m, A.m); row++)
	{
		for(int col = 0; col < fmin(n, A.n); col++)
		{
			B = mat_set(B, row, col, mat_get(A, row, col));
		}
	}
	
	return B;
}

mat_t mat_add_row(mat_t A, int i, mat_t R)
{
	if(i < 0 || i > A.m)
	{
		puts("[mat_add_row] index out of bounds");
		exit(EXIT_FAILURE);
	}
	if(R.n > A.n)
	{
		puts("[mat_add_row] incorrect operand dimensions");
		exit(EXIT_FAILURE);
	}
	if((A.m+1) * A.n > MAT_MAX_SIZE)
	{
		printf("[mat_add_row] product of matrix dimensions may not exceed %d\n", MAT_MAX_SIZE);
		exit(EXIT_FAILURE);
	}
	
	mat_t B = mat_resize(A, A.m+1, A.n);

	for(int k = B.m-2; k >= i; k--)
	{
		B = mat_set_row(B, k+1, mat_get_row(A, k)); 
	}
	B = mat_set_row(B, i, R);

	return B;
}

mat_t mat_add_col(mat_t A, int j, mat_t C)
{
	if(j < 0 || j > A.n)
	{
		puts("[mat_add_col] index out of bounds");
		exit(EXIT_FAILURE);
	}
	if(C.m > A.m)
	{
		puts("[mat_add_col] incorrect operand dimensions");
		exit(EXIT_FAILURE);
	}
	if(A.m * (A.n+1) > MAT_MAX_SIZE)
	{
		printf("[mat_add_col] product of matrix dimensions may not exceed %d\n", MAT_MAX_SIZE);
		exit(EXIT_FAILURE);
	}
	
	mat_t B = mat_resize(A, A.m, A.n+1);

	for(int k = B.m-2; k >= j; k--)
	{
		B = mat_set_col(B, k+1, mat_get_col(A, k)); 
	}
	B = mat_set_col(B, j, C);

	return B;
}

mat_t mat_del_row(mat_t A, int i)
{
	if(i < 0 || i >= A.m)
	{
		puts("[mat_del_row] index out of bounds");
		exit(EXIT_FAILURE);
	}

	mat_t B = mat_init(A.m-1, A.n);

	int j = 0;
	for(int k = 0; k < A.m; k++)
	{
		if(k != i)
		{
			B = mat_set_row(B, j, mat_get_row(A, k));
			j++;
		}
	}

	return B;
}

mat_t mat_del_col(mat_t A, int j)
{
	if(j < 0 || j >= A.n)
	{
		puts("[mat_del_col] index out of bounds");
		exit(EXIT_FAILURE);
	}

	mat_t B = mat_init(A.m, A.n-1);

	int i = 0;
	for(int k = 0; k < A.n; k++)
	{
		if(k != j)
		{
			B = mat_set_col(B, i, mat_get_col(A, k));
			i++;
		}
	}

	return B;
}

mat_t mat_swap_rows(mat_t A, int i, int j)
{
	if(i < 0 || j < 0 || i >= A.m || j >= A.m)
	{
		puts("[mat_swap_rows] index out of bounds");
		exit(EXIT_FAILURE);
	}

	mat_t temp = mat_get_row(A, i);
	A = mat_set_row(A, i, mat_get_row(A, j));
	A = mat_set_row(A, j, temp);
	return A;
}

mat_t mat_swap_cols(mat_t A, int i, int j)
{
	if(i < 0 || j < 0 || i >= A.n || j >= A.n)
	{
		puts("[mat_swap_cols] index out of bounds");
		exit(EXIT_FAILURE);
	}

	mat_t temp = mat_get_col(A, i);
	A = mat_set_col(A, i, mat_get_col(A, j));
	A = mat_set_col(A, j, temp);
	return A;
}

mat_t mat_scale(mat_t A, double l)
{
	for(int i = 0; i < A.m*A.n; i++)
	{ A.data[i] *= l; }

	return A;
}

mat_t mat_add(mat_t A, mat_t B)
{
	if(A.m != B.m || A.n != B.n)
	{
		puts("[mat_add] operand dimensions do not match");
		exit(EXIT_FAILURE);
	}

	for(int i = 0; i < A.m*A.n; i++)
	{ A.data[i] += B.data[i]; }

	return A;
}

mat_t mat_sub(mat_t A, mat_t B)
{
	if(A.m != B.m || A.n != B.n)
	{
		puts("[mat_sub] operand dimensions do not match");
		exit(EXIT_FAILURE);
	}

	for(int i = 0; i < A.m*A.n; i++)
	{ A.data[i] -= B.data[i]; }

	return A;
}

mat_t mat_mul(mat_t A, mat_t B)
{
	if(A.n != B.m)
	{
		puts("[mat_mul] incorrect operand dimensions");
		exit(EXIT_FAILURE);
	}

	mat_t C = mat_init(A.m, B.n);

	for(int i = 0; i < C.m; i++)
	{
		for(int j = 0; j < C.n; j++)
		{
			double sum = 0.0;
			for(int k = 0; k < A.n; k++)
			{
				sum += A.data[i * A.n + k] * B.data[k * B.n + j];
			}
			C.data[i * C.n + j] = sum;
			// printf("%d,%d <- %f\n", i, j, sum);
		}
	}

	return C;
}

double mat_dot(mat_t A, mat_t B)
{
	if((A.m != 1 && A.n != 1) || (B.m != 1 && B.n != 1) || (A.m * A.n != B.m * B.n))
	{
		puts("[mat_dot] operands must be vectors of equal length");
		exit(EXIT_FAILURE);
	}
	
	double sum = 0.0;
	for(int i = 0; i < A.m * A.n; i++)
	{ sum += A.data[i] * B.data[i]; }

	return sum;
}

mat_t mat_lerp(mat_t A, mat_t B, double t)
{
	if(A.m != B.m || A.n != B.n)
	{
		puts("[mat_lerp] operand dimensions do not match");
		exit(EXIT_FAILURE);
	}

	return mat_add(A, mat_scale(mat_sub(B, A), t));
}

mat_t mat_compwise(mat_t A, mat_t B)
{
	if(A.m != B.m || A.n != B.n)
	{
		puts("[mat_compwise] operand dimensions do not match");
		exit(EXIT_FAILURE);
	}

	for(int i = 0; i < A.m * A.n; i++)
	{
		B.data[i] *= A.data[i];
	}
	
	return B;
}

mat_t mat_lincomb(mat_t A, mat_t B)
{
	if(A.m != B.n)
	{
		puts("[mat_lincomb] number of first operand's rows not equal to number of second operand's columns");
		exit(EXIT_FAILURE);
	}

	mat_t S = mat_init(1, A.n);
	
	for(int i = 0; i < A.m; i++)
	{
		mat_t R = mat_scale(mat_get_row(A, i), B.data[i]);
		S = mat_add(S, R);
	}

	return S;
}

mat_t mat_trans(mat_t A)
{
	mat_t AT = mat_init(A.n, A.m);

	for(int i = 0; i < A.n; i++)
	{
		for(int j = 0; j < A.m; j++)
		{
			AT.data[i * A.m + j] = A.data[j * A.n + i];
		}
	}

	return AT;
}

double mat_norm(mat_t A)
{
	double sum = 0;
	for(int i = 0; i < A.m*A.n; i++)
	{
		sum += A.data[i] * A.data[i]; 
	}

	return sqrt(sum);
}

mat_t mat_normalized(mat_t A)
{
	return mat_scale(A, 1.0 / mat_norm(A));
}

double mat2_det(mat_t A)
{
	if(A.m != 2 || A.n != 2)
	{
		puts("[mat2_det] operand is not 2x2");
		exit(EXIT_FAILURE);
	}

	return A.data[0] * A.data[3] - A.data[1] * A.data[2]; 
}

mat_t vec3_cross(mat_t A, mat_t B)
{
	if(A.m * A.n < 3 || B.m * B.n < 3)
	{
		puts("[vec3_cross] operand is not 1x3");
		exit(EXIT_FAILURE);
	}

	double Ax = A.data[0]; double Bx = B.data[0];
	double Ay = A.data[1]; double By = B.data[1];
	double Az = A.data[2]; double Bz = B.data[2];

	return mat_row(3, Ay*Bz - Az*By, Az*Bx - Ax*Bz, Ax*By - Ay*Bx);
}

mat_t mat4_inv(mat_t A)
{
	if(A.m != 4 || A.n != 4)
	{
		puts("[mat4_inv] operand is not 4x4");
		exit(EXIT_FAILURE);
	}

	mat_t abcd = mat_del_row(A, 3);
	mat_t xyzw = mat_get_row(A, 3);
	mat_t a = mat_trans(mat_get_col(abcd, 0));
	mat_t b = mat_trans(mat_get_col(abcd, 1));
	mat_t c = mat_trans(mat_get_col(abcd, 2));
	mat_t d = mat_trans(mat_get_col(abcd, 3));
	double x = xyzw.data[0];
	double y = xyzw.data[1];
	double z = xyzw.data[2];
	double w = xyzw.data[3];
	
	mat_t s = vec3_cross(a, b);
	mat_t t = vec3_cross(c, d);
	mat_t u = mat_sub(mat_scale(a, y), mat_scale(b, x));
	mat_t v = mat_sub(mat_scale(c, w), mat_scale(d, z));

	double det = mat_dot(s, v) + mat_dot(t, u);
	mat_t row_0 = mat_add(vec3_cross(b, v), mat_scale(t, y));
	mat_t row_1 = mat_sub(vec3_cross(v, a), mat_scale(t, x));
	mat_t row_2 = mat_add(vec3_cross(d, u), mat_scale(s, w));
	mat_t row_3 = mat_sub(vec3_cross(u, c), mat_scale(s, z));
	double c00 = mat_dot(mat_scale(b, -1.0), t);
	double c01 = mat_dot(a, t);
	double c02 = mat_dot(mat_scale(d, -1.0), s);
	double c03 = mat_dot(c, s);
	mat_t col_3 = mat_col(4, c00, c01, c02, c03);

	mat_t Ainv = mat_add_row(row_0, 1, row_1);
	Ainv = mat_add_row(Ainv, 2, row_2);
	Ainv = mat_add_row(Ainv, 3, row_3);
	Ainv = mat_add_col(Ainv, 3, col_3);
	Ainv = mat_scale(Ainv, 1.0 / det);

	return Ainv;
}

mat_t scaling(double x, double y, double z)
{
	mat_t S = mat_id(4);
	S = mat_set(S, 0, 0, x);
	S = mat_set(S, 1, 1, y);
	S = mat_set(S, 2, 2, z);

	return S;
}

mat_t rotation(double x, double y, double z)
{
	mat_t Rx = mat_id(4);
	Rx = mat_set(Rx, 1, 1, cos(x));
	Rx = mat_set(Rx, 1, 2, -sin(x));
	Rx = mat_set(Rx, 2, 1, sin(x));
	Rx = mat_set(Rx, 2, 2, cos(x));

	mat_t Ry = mat_id(4);
	Ry = mat_set(Ry, 0, 0, cos(y));
	Ry = mat_set(Ry, 0, 2, sin(y));
	Ry = mat_set(Ry, 2, 0, -sin(y));
	Ry = mat_set(Ry, 2, 2, cos(y));
	mat_print(Ry);
	
	mat_t Rz = mat_id(4);
	Rz = mat_set(Rz, 0, 0, cos(z));
	Rz = mat_set(Rz, 0, 1, -sin(z));
	Rz = mat_set(Rz, 1, 0, sin(z));
	Rz = mat_set(Rz, 1, 1, cos(z));

	return mat_mul(Rz, mat_mul(Ry, Rx));
}

mat_t translation(double x, double y, double z)
{
	mat_t T = mat_id(4);
	mat_t col = mat_col(4, x, y, z, 1.0);
	T = mat_set_col(T, 3, col);

	return T;
}

mat_t lookat(mat_t eye, mat_t target, mat_t up)
{
	eye = mat_scale(mat_resize(eye, 1, 3), -1.0);
	target = mat_resize(target, 1, 3);
	up = mat_normalized(mat_resize(up, 1, 3));

	mat_t w = mat_normalized(mat_sub(target, eye));
	mat_t u = mat_normalized(vec3_cross(up, w));
	mat_t v = mat_normalized(vec3_cross(w, u));

	mat_t A = mat_id(4);
	A = mat_blit(A, 0, 0, u);
	A = mat_blit(A, 1, 0, v);
	A = mat_blit(A, 2, 0, w);
	
	mat_t T = translation(eye.data[0], eye.data[1], eye.data[2]);

	return mat_mul(A, T);
}

mat_t perspective(double aspect, double fovy, double near, double far)
{
	double foc_len = 1.0 / tan(fovy * 0.5);
	double z_norm = (far + near) / (far - near);
	double z_off = (-2.0 * far * near) / (far - near);

	mat_t P = mat_init(4,4);
	P = mat_set(P, 0, 0,  foc_len / aspect);
	P = mat_set(P, 1, 1, foc_len);
	P = mat_set(P, 2, 2, z_norm);
	P = mat_set(P, 2, 3, z_off);
	P = mat_set(P, 3, 2, -1.0);

	return P;
}

mat_t window(int width, int height)
{
	mat_t S = mat_id(4);
	S = mat_set(S, 0, 0, width * 0.5);
	S = mat_set(S, 0, 3, width * 0.5);
	S = mat_set(S, 1, 1, height * 0.5);
	S = mat_set(S, 1, 3, height * 0.5);
	S = mat_set(S, 2, 2, 0.5);
	S = mat_set(S, 2, 3, 0.5);

	return S;
}

void mat_print(mat_t A)
{
	printf("Matrix %d x %d\n", A.m, A.n);
	for(int i = 0; i < A.m; i++)
	{
		printf("[");
		for(int j = 0; j < A.n; j++)
		{
			printf(" %.6e", A.data[i * A.n + j]);
		}
		printf(" ]\n");
	}
	printf("\n");
}

void mat_dump(mat_t A)
{
	printf("[RAW] Matrix %d x %d\n", A.m, A.n);
	printf("[ ");
	for(int i = 0; i < MAT_MAX_SIZE; i++)
	{
		printf("%.6e", A.data[i]);
		if(i < MAT_MAX_SIZE-1)
		{
			printf(", ");
		}
	}
	printf(" ]\n");
}

void num_print(double n)
{
	printf("%10.10f\n", n);
}
