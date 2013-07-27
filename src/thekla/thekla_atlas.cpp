
#include "thekla_atlas.h"

#include "nvmesh/halfedge/Edge.h"
#include "nvmesh/halfedge/Mesh.h"
#include "nvmesh/halfedge/Face.h"
#include "nvmesh/halfedge/Vertex.h"
#include "nvmesh/param/Atlas.h"

#include "nvmath/Vector.inl"

#include "nvcore/Array.inl"


using namespace Thekla;
using namespace nv;


static Atlas_Error input_to_mesh(const Atlas_Input_Mesh * input, HalfEdge::Mesh * mesh) {

    Array<uint> canonicalMap;
    canonicalMap.reserve(input->vertex_count);

    for (int i = 0; i < input->vertex_count; i++) {
        const Atlas_Input_Vertex & input_vertex = input->vertex_array[i];
        const float * pos = input_vertex.position;
        const float * nor = input_vertex.normal;
        const float * tex = input_vertex.uv;

        HalfEdge::Vertex * vertex = mesh->addVertex(Vector3(pos[0], pos[1], pos[2]));
        vertex->nor.set(nor[0], nor[1], nor[2]);
        vertex->tex.set(tex[0], tex[1]);

        canonicalMap.append(input_vertex.first_colocal);
    }

    mesh->linkColocalsWithCanonicalMap(canonicalMap);


    const int face_count = input->face_count;

    int non_manifold_faces = 0;
    for (int i = 0; i < face_count; i++) {
        const Atlas_Input_Face & input_face = input->face_array[i];

        int v0 = input_face.vertex_index[0];
        int v1 = input_face.vertex_index[1];
        int v2 = input_face.vertex_index[2];

        HalfEdge::Face * face = mesh->addFace(v0, v1, v2);
        if (face != NULL) {
            face->material = input_face.material_index;
        }
        else {
            non_manifold_faces++;
        }
    }

    mesh->linkBoundary();

    /*if (non_manifold_faces != 0) {
        return Atlas_Error_Invalid_Mesh_Non_Manifold;
    }*/

    return Atlas_Error_Success;
}

static Atlas_Error mesh_atlas_to_output(const HalfEdge::Mesh * mesh, const Atlas & atlas, Atlas_Output_Mesh * output) {

    /*
    struct Atlas_Output_Vertex {
        float uv[2];
        int xref;       // Index of vertex from which this vertex originated.
    };

    struct Atlas_Output_Mesh {
        int vertex_count;
        int index_count;
        Atlas_Output_Vertex * vertex_array;
        int * index_array;
    };
    */

    // Allocate vertices.
    const int vertex_count = atlas.vertexCount();
    output->vertex_count = vertex_count;
    output->vertex_array = new Atlas_Output_Vertex[vertex_count];

    // Output vertices.
    const int chart_count = atlas.chartCount();
    for (int i = 0; i < chart_count; i++) {
        const Chart * chart = atlas.chartAt(i);
        uint vertexOffset = atlas.vertexCountBeforeChartAt(i);

        const uint chart_vertex_count = chart->vertexCount();
        for (uint v = 0; v < chart_vertex_count; v++) {
            Atlas_Output_Vertex & output_vertex = output->vertex_array[vertexOffset + v]; 

            uint original_vertex = chart->mapChartVertexToOriginalVertex(v);
            output_vertex.xref = original_vertex;

            Vector2 uv = chart->chartMesh()->vertexAt(v)->tex;
            output_vertex.uv[0] = uv.x;
            output_vertex.uv[1] = uv.y;
        }
    }

    const int face_count = mesh->faceCount();
    output->index_count = face_count * 3;
    output->index_array = new int[face_count * 3];

        // Set face indices.
    for (int f = 0; f < face_count; f++) {
        uint c = atlas.faceChartAt(f);
        uint i = atlas.faceIndexWithinChartAt(f);
        uint vertexOffset = atlas.vertexCountBeforeChartAt(c);

        const Chart * chart = atlas.chartAt(c);
        nvDebugCheck(chart->faceAt(i) == f);

        const HalfEdge::Face * face = chart->chartMesh()->faceAt(i);
        const HalfEdge::Edge * edge = face->edge;

        output->index_array[3*f+0] = vertexOffset + edge->vertex->id;
        output->index_array[3*f+1] = vertexOffset + edge->next->vertex->id;
        output->index_array[3*f+2] = vertexOffset + edge->next->next->vertex->id;
    }

    return Atlas_Error_Success;
}



extern "C" Atlas_Error atlas_generate(const Atlas_Options * options, const Atlas_Input_Mesh * input, Atlas_Output_Mesh * output) {
    // Validate args.
    if (options == NULL || input == NULL || output == NULL) return Atlas_Error_Invalid_Args;

    // Validate options.
    if (options->charter != Atlas_Charter_Extract && options->charter != Atlas_Charter_Witness) {
        return Atlas_Error_Invalid_Options;
    }
    if (options->charter == Atlas_Charter_Extract && options->charter_option_count != 0) {
        return Atlas_Error_Invalid_Options;
    }
    if (options->charter == Atlas_Charter_Witness && options->charter_option_count != 1) {
        return Atlas_Error_Invalid_Options;
    }

    if (options->mapper != Atlas_Mapper_LSCM) {
        return Atlas_Error_Invalid_Options;
    }
    if (options->mapper == Atlas_Mapper_LSCM && options->charter_option_count != 0) {
        return Atlas_Error_Invalid_Options;
    }

    if (options->packer != Atlas_Packer_Witness) {
        return Atlas_Error_Invalid_Options;
    }
    if (options->packer == Atlas_Packer_Witness && options->charter_option_count != 1) {
        return Atlas_Error_Invalid_Options;
    }

    // Validate input mesh.
    for (int i = 0; i < input->face_count; i++) {
        int v0 = input->face_array[i].vertex_index[0];
        int v1 = input->face_array[i].vertex_index[1];
        int v2 = input->face_array[i].vertex_index[2];

        if (v0 < 0 || v0 >= input->vertex_count || 
            v1 < 0 || v1 >= input->vertex_count || 
            v2 < 0 || v2 >= input->vertex_count)
        {
            return Atlas_Error_Invalid_Mesh;
        }
    }


    // Build half edge mesh.
    AutoPtr<HalfEdge::Mesh> mesh(new HalfEdge::Mesh);

    Atlas_Error error = input_to_mesh(input, mesh.ptr());
    if (error != Atlas_Error_Success) {
        return error;
    }


    // Charter.
    Atlas atlas(mesh.ptr());
    
    if (options->charter == Atlas_Charter_Extract) {
        return Atlas_Error_Not_Implemented;
    }
    else if (options->charter == Atlas_Charter_Witness) {
        // @@ Init segmentations settings!
        SegmentationSettings segmentation_settings;
        // options->charter_options[0];

        atlas.computeCharts(segmentation_settings);
    }


    // Mapper.
    if (options->mapper == Atlas_Mapper_LSCM) {
        atlas.parameterizeCharts();
    }


    // Packer.
    if (options->packer == Atlas_Packer_Witness) {
        int packing_quality = 0; // best
        float texel_area = options->packer_options[0];
        int texel_padding = 1;      // @@ Add more padding and block alignment options.

        /*float utilization =*/ atlas.packCharts(packing_quality, texel_area, texel_padding);
    }


    // Build output mesh.
    return mesh_atlas_to_output(mesh.ptr(), atlas, output);
}
