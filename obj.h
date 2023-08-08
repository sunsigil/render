#ifndef OBJ_H
#define OBJ_H

#include <stdio.h>
#include "linalg.h"

typedef struct OBJ
{
	size_t file_size;

	double* vs;
	double* vts;
	double* vns;
	int* fs;

	int v_count;
	int vt_count;
	int vn_count;
	int f_count;

} OBJ_t;

OBJ_t OBJ_read(char* path);
void OBJ_extract_vs(OBJ_t obj, mat_t* vs);
void OBJ_extract_vts(OBJ_t obj, mat_t* vts);
void OBJ_extract_vns(OBJ_t obj, mat_t* vns);
void OBJ_extract_fs(OBJ_t obj, mat_t* fs);
mat_t OBJ_face_vs(OBJ_t obj, int f_idx);
mat_t OBJ_face_vts(OBJ_t obj, int f_idx);
mat_t OBJ_face_vns(OBJ_t obj, int f_idx);
void OBJ_dispose(OBJ_t* obj);

#endif
