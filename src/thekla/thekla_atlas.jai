
// Thekla Atlas Generator

Thekla :: struct {
    Atlas_Charter :: enum u32 {
        Witness;
        Extract;
        Default :: Witness;
    }

    Atlas_Mapper :: enum u32 {
        LSCM;
        Default :: LSCM;
    }

    Atlas_Packer :: enum u32 {
        Witness;
        Default :: Witness;
    }

    Atlas_Options :: struct {
        charter: Atlas_Charter;
        charter_options: union {
            witness: struct {
                proxy_fit_metric_weight: float;
                roundness_metric_weight: float;
                straightness_metric_weight: float;
                normal_seam_metric_weight: float;
                texture_seam_metric_weight: float;
                max_chart_area: float;
                max_boundary_length: float;
            };
            extract: struct {
            };
        };

        mapper: Atlas_Mapper;
        mapper_options: struct {
            preserve_uvs: bool;
            preserve_boundary: bool;
        };

        packer: Atlas_Packer;
        packer_options: union {
            witness: struct {
                packing_quality: s32;
                texels_per_unit: float; // Unit to texel scale. e.g. a 1x1 quad with texelsPerUnit of 32 will take up approximately 32x32 texels in the atlas.
                block_align: bool;      // Align charts to 4x4 blocks. 
                conservative: bool;     // Pack charts with extra padding.
            };
        };
    }
    #assert(size_of(Atlas_Options) == 56);


    Atlas_Input_Vertex :: struct {
        position: [3]float;
        normal: [3]float;
        uv: [2]float;
        first_colocal: s32;
    }

    Atlas_Input_Face :: struct {
        vertex_index: [3]s32;
        material_index: s32;
    }

    Atlas_Input_Mesh :: struct {
        vertex_count: s32;
        face_count: s32;
        vertex_array: *Atlas_Input_Vertex;
        face_array: *Atlas_Input_Face;
    }

    Atlas_Output_Vertex :: struct {
        uv: [2]float;
        xref: s32;   // Index of input vertex from which this output vertex originated.
    }

    Atlas_Output_Mesh :: struct {
        atlas_width: s32;
        atlas_height: s32;
        vertex_count: s32;
        index_count: s32;
        vertex_array: *Atlas_Output_Vertex;
        index_array: *s32;
    }

    Atlas_Error :: enum {
        Success;
        Invalid_Args;
        Invalid_Options;
        Invalid_Mesh;
        Invalid_Mesh_Non_Manifold;
        Not_Implemented;
    }

    atlas_set_default_options :: (options: *Atlas_Options) #foreign thekla_atlas "?atlas_set_default_options@Thekla@@YAXPEAUAtlas_Options@1@@Z";
    atlas_generate :: (input: *Atlas_Input_Mesh, options: *Atlas_Options, error: *Atlas_Error) -> *Atlas_Output_Mesh #foreign thekla_atlas "?atlas_generate@Thekla@@YAPEAUAtlas_Output_Mesh@1@PEBUAtlas_Input_Mesh@1@PEBUAtlas_Options@1@PEAW4Atlas_Error@1@@Z";
    atlas_free :: (output: *Atlas_Output_Mesh) #foreign thekla_atlas "?atlas_free@Thekla@@YAXPEAUAtlas_Output_Mesh@1@@Z";

} // Thekla;


#scope_file 

thekla_atlas :: #foreign_library "thekla_atlas";



