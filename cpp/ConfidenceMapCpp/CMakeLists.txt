CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )
#all headers
SET(headers ConfidenceMaps2D.h
			ConfidenceMaps2DFacade.h
			RandomWalks2D.h
			RandomWalks2DFacade.h
			RandomWalks3D.h
			RandomWalks3DFacade.h
			RandomWalksCore.h
			SparseSolverEigenCG.h
			SparseSolverEigenCustom.h
			SparseSolverEigenLLT.h
			SparseSolverFactory.h
			SparseSolverInterface.h
			SparseSolverViennaCPU.h
			SparseSolverViennaGPU.h)

SET(sources ConfidenceMaps2D.cpp
			ConfidenceMaps2DFacade.cpp
			RandomWalks2D.cpp
			RandomWalks2DFacade.cpp
			RandomWalks3D.cpp
			RandomWalks3DFacade.cpp
			RandomWalksCore.cpp
			SparseSolverEigenCG.cpp
			SparseSolverEigenCustom.cpp
			SparseSolverEigenLLT.cpp
			SparseSolverFactory.cpp
			SparseSolverViennaCPU.cpp
			SparseSolverViennaGPU.cpp)
#include eigen and ViennaCl headers
INCLUDE_DIRECTORIES(extra_part/eigen3)
INCLUDE_DIRECTORIES(extra_part/ViennaCL-1.7.1)
INCLUDE_DIRECTORIES(C:/Program\ Files/NVIDIA\ GPU\ Computing\ Toolkit/CUDA/v8.0/include)
INCLUDE_DIRECTORIES(C:/Program\ Files/NVIDIA\ GPU\ Computing\ Toolkit/CUDA/v8.0/include/CL)
#path to OpenCl binaries
LINK_DIRECTORIES(C:/Program\ Files/NVIDIA\ GPU\ Computing\ Toolkit/CUDA/v8.0/lib/x64)
ADD_LIBRARY( RandomWalksLib STATIC ${sources} ${headers} )

##put debug target to x64/Debug and all other configurations to x64/Release
#SET_TARGET_PROPERTIES( RandomWalksLib PROPERTIES
#						ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/x64/Release #static libs are archives
#						ARCHIVE_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/x64/Debug
#						ARCHIVE_OUTPUT_NAME RandomWalksLib
#						ARCHIVE_OUTPUT_NAME_DEBUG RandomWalksLibd ) #debug build should end with a 'd'
TARGET_LINK_LIBRARIES( RandomWalksLib OpenCL )
