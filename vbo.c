#include "vbo.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>

void vert_node_dispose(vert_node_t* node)
{
	if(node == NULL)
	{ return; }
	vert_node_dispose(node->next);
	free(node);
}

vert_map_t vert_map_init(int n)
{
	vert_map_t map;
	map.verts = calloc(n, sizeof(vert_node_t*));
	map.capacity = n;
	map.count = 0;
	return map;
}

unsigned int hash(int v, int vt, int vn, int n)
{
	unsigned int h = 1;
	h = 31 * h + v;
	h = 31 * h + vt;
	h = 31 * h + vn;
	return h % n;
}

int vert_map_add(vert_map_t* map, int v, int vt, int vn)
{
	vert_node_t* vert = malloc(sizeof(vert_node_t));
	vert->v = v;
	vert->vt = vt;
	vert->vn = vn;
	vert->idx = map->count;
	vert->next = NULL;

	int idx = hash(v, vt, vn, map->capacity);
	if(map->verts[idx] != NULL)
	{ vert->next = map->verts[idx]; }

	map->verts[idx] = vert;
	map->count++;

	return vert->idx;
}

int vert_map_lookup(vert_map_t* map, int v, int vt, int vn)
{
	int idx = hash(v, vt, vn, map->capacity);
	vert_node_t* ptr = map->verts[idx];
	while(ptr != NULL)
	{
		if(ptr->v == v && ptr->vt == vt && ptr->vn == vn)
		{ return ptr->idx; }
		ptr = ptr->next;
	}
	return -1;
}

void vert_map_dispose(vert_map_t map)
{
	for(int i = 0; i < map.capacity; i++)
	{ vert_node_dispose(map.verts[i]); }
	free(map.verts);
}

VBO_t VBO_init(OBJ_t obj)
{
	vert_map_t map = vert_map_init(obj.f_count / 4);
	
	int* indices = malloc(sizeof(int) * obj.f_count*3);
	for(int i = 0; i < obj.f_count; i++)
	{
		int fidx = i*9;
		int iidx = i*3;
		
		for(int j = 0; j < 3; j++)
		{
			int pidx = fidx+j*3;

			int v = obj.fs[pidx+0];
			int vt = obj.fs[pidx+1];
			int vn = obj.fs[pidx+2];

			int vidx = vert_map_lookup(&map, v, vt, vn);
			if(vidx == -1)
			{ vidx = vert_map_add(&map, v, vt, vn); }

			indices[iidx+j] = vidx;
		}
	}
	
	double* vertices = malloc(sizeof(double) * map.count*8);
	for(int i = 0; i < map.capacity; i++)
	{
		vert_node_t* ptr = map.verts[i];

		while(ptr != NULL)
		{
			int vidx = ptr->idx*8;

			int v = ptr->v * 3;
			int vt = ptr->vt * 2;
			int vn = ptr->vn * 3;

			vertices[vidx+0] = obj.vs[v+0];
			vertices[vidx+1] = obj.vs[v+1];
			vertices[vidx+2] = obj.vs[v+2];

			vertices[vidx+3] = obj.vts[vt+0];
			vertices[vidx+4] = obj.vts[vt+1];
			
			vertices[vidx+5] = obj.vns[vn+0];
			vertices[vidx+6] = obj.vns[vn+1];
			vertices[vidx+7] = obj.vns[vn+2];

			ptr = ptr->next;
		}
	}

	vert_map_dispose(map);

	VBO_t vbo;
	vbo.vertices = vertices;
	vbo.indices = indices;
	vbo.face_count = obj.f_count;
	
	return vbo;
}

void VBO_dispose(VBO_t vbo)
{
	free(vbo.vertices);
	free(vbo.indices);
}
