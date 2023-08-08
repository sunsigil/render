#include "obj.h"

typedef struct vert_node
{
	int v;
	int vt;
	int vn;
	int idx;
	struct vert_node* next;
} vert_node_t;

typedef struct vert_map
{
	vert_node_t** verts;
	int capacity;
	int count;
} vert_map_t;

typedef struct VBO
{
	double* vertices;
	int* indices;
	int face_count;
} VBO_t;

void vert_node_dispose(vert_node_t* node);
vert_map_t vert_map_init(int n);
int vert_map_add(vert_map_t* map, int v, int vt, int vn);
int vert_map_lookup(vert_map_t* map, int v, int vt, int vn);
void vert_map_dispose(vert_map_t map);
VBO_t VBO_init(OBJ_t obj);
void VBO_dispose(VBO_t vbo);

