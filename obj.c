#include "obj.h"

#include "linalg.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>
#include <math.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>

int is_obj_num(char c)
{
	return
	isdigit(c) ||
	c == '.' ||
	c == '-' ||
	c == 'e';
}

void parse_line(char* line, double* dest)
{
	int pos = 0;
	int vals = 0;

	while(line[pos] != '\n' && line[pos] != EOF)
	{
		if(is_obj_num(line[pos]))
		{
			double d = atof(line+pos);
			dest[vals] = d;
			vals++;
			
			while(is_obj_num(line[pos]))
			{ pos++; }
		}
		else
		{ pos++; }
	}
}

OBJ_t OBJ_read(char* path)
{
	int fd = open(path, O_RDONLY, S_IRUSR);
	if(fd == -1)
	{
		perror("[OBJ_read] open");
		exit(EXIT_FAILURE);
	}
	
	uint32_t byte_count = lseek(fd, 0, SEEK_END);
	lseek(fd, 0, SEEK_SET);
	uint8_t* bytes = mmap(NULL, byte_count, PROT_READ, MAP_PRIVATE, fd, SEEK_SET);
	if(bytes == NULL)
	{
		perror("[OBJ_read] mmap");
		exit(EXIT_FAILURE);
	}
	
	int v_count = 0;
	int vt_count = 0;
	int vn_count = 0;
	int f_count = 0;

	int pos = 0;
	while(pos < byte_count)
	{
		if(memcmp(bytes+pos, "v ", 2) == 0)
		{ v_count++; }
		else if(memcmp(bytes+pos, "vt ", 3) == 0)
		{ vt_count++; }
		else if(memcmp(bytes+pos, "vn ", 3) == 0)
		{ vn_count++; }
		else if(memcmp(bytes+pos, "f ", 2) == 0)
		{ f_count++; }
		
		while(pos < byte_count && bytes[pos] != '\n')
		{ pos++; }
		pos++;
	}

	double* vs = malloc(v_count*3 * sizeof(double));
	int vs_read = 0;
	double* vts = malloc(vt_count*2 * sizeof(double));
	int vts_read = 0;
	double* vns = malloc(vn_count*3 * sizeof(double));
	int vns_read = 0;
	double* fs = malloc(f_count*9 * sizeof(double));
	int fs_read = 0;

	pos = 0;
	while(pos < byte_count)
	{
		if(memcmp(bytes+pos, "v ", 2) == 0)
		{
			parse_line(bytes+pos, vs+vs_read);
			vs_read += 3;
		}
		else if(memcmp(bytes+pos, "vt ", 3) == 0)
		{
			parse_line(bytes+pos, vts+vts_read);
			vts_read += 2;
		}
		else if(memcmp(bytes+pos, "vn ", 3) == 0)
		{
			parse_line(bytes+pos, vns+vns_read);
			vns_read += 3;
		}
		else if(memcmp(bytes+pos, "f ", 2) == 0)
		{
			parse_line(bytes+pos, fs+fs_read);
			fs_read += 9;
		}
		
		while(pos < byte_count && bytes[pos] != '\n')
		{ pos++; }
		pos++;
	}
	
	OBJ_t obj;
	obj.file_size = byte_count;

	obj.vs = vs;
	obj.vts = vts;
	obj.vns = vns;

	obj.fs = malloc(f_count*9 * sizeof(int));
	for(int i = 0; i < f_count * 9; i++)
	{ obj.fs[i] = ((int) fs[i]) - 1; }
	free(fs);

	obj.v_count = v_count;
	obj.vt_count = vt_count;
	obj.vn_count = vn_count;
	obj.f_count = f_count;
	
	if(munmap(bytes, byte_count) == -1)
	{
		perror("[OBJ_read] munmap");
		exit(EXIT_FAILURE);
	}
	if(close(fd) == -1)
	{
		perror("[OBJ_read] close");
		exit(EXIT_FAILURE);
	}

	return obj;
}

void OBJ_extract_vs(OBJ_t obj, mat_t* vs)
{
	for(int i = 0; i < obj.v_count*3; i+=3)
	{ vs[i/3] = mat_row(4, obj.vs[i], obj.vs[i+1], obj.vs[i+2], 1.0); }
}

void OBJ_extract_vts(OBJ_t obj, mat_t* vts)
{
	for(int i = 0; i < obj.vt_count*2; i+=2)
	{ vts[i/2] = mat_row(2, obj.vts[i], obj.vts[i+1]); }
}

void OBJ_extract_vns(OBJ_t obj, mat_t* vns)
{
	for(int i = 0; i < obj.vn_count*3; i+=3)
	{ vns[i/3] = mat_row(4, obj.vns[i], obj.vns[i+1], obj.vns[i+2], 0.0); }
}

void OBJ_extract_fs(OBJ_t obj, mat_t* fs)
{
	for(int i = 0; i < obj.f_count*9; i+=9)
	{
		mat_t a = mat_row(3, (double) obj.fs[i+0], (double) obj.fs[i+1], (double) obj.fs[i+2]);
		mat_t b = mat_row(3, (double) obj.fs[i+3], (double) obj.fs[i+4], (double) obj.fs[i+5]);
		mat_t c = mat_row(3, (double) obj.fs[i+6], (double) obj.fs[i+7], (double) obj.fs[i+8]);
		mat_t f = mat_init(3, 3);
		f = mat_set_row(f, 0, a);
		f = mat_set_row(f, 1, b);
		f = mat_set_row(f, 2, c);
		fs[i/9] = f;
	}
}

mat_t OBJ_face_vs(OBJ_t obj, int f_idx)
{
	int* face = obj.fs + (f_idx*9);
	mat_t vs = mat_init(3, 4);

	for(int i = 0; i < 3; i++)
	{
		int v_idx = face[i*3+0] * 3;
		mat_t v = mat_row(4, obj.vs[v_idx+0], obj.vs[v_idx+1], obj.vs[v_idx+2], 1.0);
		vs = mat_set_row(vs, i, v);
	}

	return vs;
}

mat_t OBJ_face_vts(OBJ_t obj, int f_idx)
{
	int* face = obj.fs + (f_idx*9);
	mat_t vts = mat_init(3, 2);

	for(int i = 0; i < 3; i++)
	{
		int vt_idx = face[i*3+1] * 2;
		mat_t vt = mat_row(2, obj.vts[vt_idx+0], obj.vts[vt_idx+1]);
		vts = mat_set_row(vts, i, vt);
	}

	return vts;
}

mat_t OBJ_face_vns(OBJ_t obj, int f_idx)
{
	int* face = obj.fs + (f_idx*9);
	mat_t vns = mat_init(3, 4);

	for(int i = 0; i < 3; i++)
	{
		int vn_idx = face[i*3+2] * 3;
		mat_t vn = mat_row(4, obj.vns[vn_idx+0], obj.vns[vn_idx+1], obj.vns[vn_idx+2], 0.0);
		vns = mat_set_row(vns, i, vn);
	}

	return vns;
}

void OBJ_dispose(OBJ_t* obj)
{
	free(obj->vs);
	free(obj->vns);
	free(obj->vts);
	free(obj->fs);
}

