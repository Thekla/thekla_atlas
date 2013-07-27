PROJECT(nvmesh)

SET(MESH_SRCS
	nvmesh.h
	BaseMesh.h
	BaseMesh.cpp
	MeshAdjacency.h
	MeshAdjacency.cpp
	MeshBuilder.h
	MeshBuilder.cpp
	TriMesh.h
	TriMesh.cpp
	QuadTriMesh.h
	QuadTriMesh.cpp
	PolyMesh.h
	VertexData.h
	VertexData.cpp
	MeshTopology.h
	MeshTopology.cpp
	MeshSmoothGroup.h
	halfedge/HalfEdge.h
	halfedge/HalfEdge.cpp
	halfedge/HalfEdgeMesh.h
	halfedge/HalfEdgeMesh.cpp
	halfedge/HalfEdgeFace.h
	halfedge/HalfEdgeFace.cpp
	halfedge/HalfEdgeVertex.h
	subdiv/LoopGregoryApproximation.h
	subdiv/LoopGregoryApproximation.cpp
	subdiv/CubicGregoryMesh.h
	subdiv/CubicGregoryMesh.cpp
	subdiv/Subdivide.h
	subdiv/Subdivide.cpp
	subdiv/RemapFaces.h
	subdiv/RemapFaces.cpp
	subdiv/AccMesh.h
	subdiv/AccMesh.cpp
	subdiv/AccPatch.h
	subdiv/AccPatch.cpp
	subdiv/AccMeshBuilder.h
	subdiv/AccMeshBuilder.cpp
	subdiv/StencilMask.h
	subdiv/StencilMask.cpp
	import/MeshImport.h
	import/MeshImport.cpp
	import/MeshImportOBJ.h
	import/MeshImportOBJ.cpp
	import/MeshImportM.h
	import/MeshImportM.cpp
	import/MeshImportMD5.h
	import/MeshImportMD5.cpp
	import/UnrealFile.h
	import/UnrealFile.cpp
	import/MeshImportPSK.h
	import/MeshImportPSK.cpp
	import/MeshImportOFF.h
	import/MeshImportOFF.cpp
	export/MeshExportOBJ.h
	export/MeshExportOBJ.cpp
	animation/MeshSkeleton.h
	animation/MeshSkeleton.cpp
	animation/Bone.h
	animation/Animation.h
	animation/Animation.cpp
	animation/CudaSkinning.h
	animation/CudaSkinning.cpp
	animation/MorphTarget.h
	animation/MorphTarget.cpp
	geometry/Bounds.h
	geometry/Bounds.cpp
	geometry/MeshTransform.h
	geometry/MeshTransform.cpp
	geometry/MeshNormals.h
	geometry/MeshNormals.cpp
	geometry/TangentSpace.h
	geometry/TangentSpace.cpp
	raster/Raster.h
	raster/Raster.cpp
	raster/ClippedTriangle.h
	raster/ClippedTriangle.cpp
	render/VertexCache.h
	render/VertexCache.cpp
	render/MeshOptimizer.h
	render/MeshOptimizer.cpp
	render/NvTriStrip.h
	render/NvTriStrip.cpp
	render/NvTriStripObjects.h
	render/NvTriStripObjects.cpp
	render/NvTriStripVertexCache.h
	render/NvTriStripVertexCache.cpp
	param/Atlas.h
	param/Atlas.cpp
	param/Seams.h
	param/Seams.cpp
	param/LeastSquaresConformalMap.h
	param/LeastSquaresConformalMap.cpp
	param/ParameterizationQuality.h
	param/ParameterizationQuality.cpp
	param/SingleFaceMap.h
	param/SingleFaceMap.cpp
	mender/NVMeshMender.h
	mender/NVMeshMender.cpp
	weld/VertexWeld.h
	weld/VertexWeld.cpp
	weld/Weld.h
	weld/Snap.h
	weld/Snap.cpp
	raytracing/Raytracing.h
	raytracing/FaceBuffer.h
	raytracing/FaceBuffer.cpp
	raytracing/Grid.h
	raytracing/Grid.cpp
	kdtree/MeshKDTree.h
	kdtree/MeshKDTree.cpp
	kdtree/KDTree.h
	kdtree/KDTree.cpp
	kdtree/KDTreeBuilder.h
	kdtree/KDTreeBuilder.cpp
	qms/QuadSimplificationMesh.h
	qms/QuadSimplificationMesh.cpp)

#SUBDIRS(import)
#SUBDIRS(geometry)
#SUBDIRS(halfedge)
#SUBDIRS(subdiv)
#SUBDIRS(viewer)
SUBDIRS(test)
SUBDIRS(tools)

#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# if Tootle
#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/render/tootle/include)
#SET(LIBS ${LIBS} ${CMAKE_CURRENT_SOURCE_DIR}/render/tootle/lib/libTootle.a)

# targets
ADD_DEFINITIONS(-DNVMESH_EXPORTS)

IF(NVMESH_SHARED)	
	ADD_LIBRARY(nvmesh SHARED ${MESH_SRCS})
ELSE(NVMESH_SHARED)
	ADD_LIBRARY(nvmesh ${MESH_SRCS})
ENDIF(NVMESH_SHARED)

TARGET_LINK_LIBRARIES(nvmesh ${LIBS} nvcore nvmath nvimage)

INSTALL(TARGETS nvmesh
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib/static)

