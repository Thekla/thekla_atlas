
// Thekla Atlas Generator

#ifdef __cplusplus
namespace Thekla {
#endif


enum Atlas_Charter {
    Atlas_Charter_Extract,  // Options: ---
    Atlas_Charter_Witness,	// Options: threshold
    Atlas_Charter_Default = Atlas_Charter_Witness
};

enum Atlas_Mapper {
    Atlas_Mapper_LSCM,      // Options: ---
    Atlas_Mapper_Default = Atlas_Mapper_LSCM
};

enum Atlas_Packer {
    Atlas_Packer_Witness,	// Options: texel_area
    Atlas_Packer_Default = Atlas_Packer_Witness
};

struct Atlas_Options {
    Atlas_Charter charter;
    float * charter_options;
   	int charter_option_count;

    Atlas_Mapper mapper;
    float * mapper_options;
    int mapper_opion_count;

    Atlas_Packer packer;
    float * packer_options;
    int packer_option_count;
};

struct Atlas_Input_Vertex {
    float position[3];
    float normal[3];
    float uv[2];
    int first_colocal;
};

struct Atlas_Input_Face {
    int vertex_index[3];
    int material_index;
};

struct Atlas_Input_Mesh {
    int vertex_count;
    int face_count;
    Atlas_Input_Vertex * vertex_array;
    Atlas_Input_Face * face_array;
};

struct Atlas_Output_Vertex {
    float uv[2];
    int xref;   // Index of input vertex from which this output vertex originated.
};

struct Atlas_Output_Mesh {
    int vertex_count;
    int index_count;
    Atlas_Output_Vertex * vertex_array;
    int * index_array;
};

enum Atlas_Error {
    Atlas_Error_Success,
    Atlas_Error_Invalid_Args,
    Atlas_Error_Invalid_Options,
    Atlas_Error_Invalid_Mesh,
    Atlas_Error_Invalid_Mesh_Non_Manifold,
    Atlas_Error_Not_Implemented,
};

extern "C" Atlas_Error atlas_generate(const Atlas_Options * options, const Atlas_Input_Mesh * input, Atlas_Output_Mesh * output);


//extern "C" Atlas_Input_Mesh * mesh_read_obj(const char * file_name);


/*

Should we represent the input mesh with an opaque structure that simply holds pointers to the user data? That would allow us to avoid having to copy attributes to an intermediate representation.

struct Atlas_Input_Mesh;

void mesh_set_vertex_position(Atlas_Input_Mesh * mesh, float * ptr, int stride);
void mesh_set_vertex_normal(Atlas_Input_Mesh * mesh, float * ptr, int stride);
void mesh_set_vertex_uv(Mesh * mesh, float * ptr, int stride);

void mesh_set_index(Mesh * mesh, int * ptr);
*/

#ifdef __cplusplus
} // thekla namespace
#endif

