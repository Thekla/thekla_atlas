
#include <stdio.h> // FILE

namespace Thekla {

struct Obj_Vertex {
	float position[3];
	float normal[3];
	float uv[2];
	int first_colocal;
};

struct Obj_Face {
	int vertex_index[3];
	int material_index;
};

struct Obj_Material {
	// @@ Read obj mtl parameters as is.
};

struct Obj_Mesh {
	int vertex_count;
	Obj_Vertex * vertex_array;

	int face_count;
	Obj_Face * face_array;

	int material_count;
	Obj_Material * material_array;
};

enum Load_Flags {
	Load_Flag_Weld_Attributes,
};

struct Obj_Load_Options {
	int load_flags;
};

Obj_Mesh * obj_mesh_load(const char * filename, const Obj_Load_Options * options);
Obj_Mesh * obj_mesh_load_from_file(FILE * file, const Obj_Load_Options * options);
void obj_mesh_free(Obj_Mesh * mesh);



Obj_Mesh * obj_mesh_load(const char * filename, const Obj_Load_Options * options) {
    Obj_Mesh * output = NULL;

    if (FILE * fp = fopen(filename, "rb")) {
        output = obj_mesh_load_from_file(fp, options);
        fclose(fp);
    }

    return output;
}

Obj_Mesh * obj_mesh_load_from_file(FILE * file, const Obj_Load_Options * options) {

    // First, read as is with minimal translation.

    // Then, transform to output representation.

    return NULL;
}



void obj_mesh_free(Obj_Mesh * mesh) {
    if (mesh != NULL) {
        delete [] mesh->vertex_array;
        delete [] mesh->face_array;
        delete [] mesh->material_array;
        delete mesh;
    }
}

} // Thekla namespace
