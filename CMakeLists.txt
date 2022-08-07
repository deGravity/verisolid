cmake_minimum_required(VERSION 3.8)
project(verisolid)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

#find_package(pybind11 CONFIG REQUIRED)
find_package(Eigen3 REQUIRED)
#find_package(OpenCASCADE REQUIRED)

add_subdirectory(parasolid)

add_definitions(-DTEST_X_T="${CMAKE_CURRENT_SOURCE_DIR}/brick.x_t")
add_definitions(-DTEST_STEP="${CMAKE_CURRENT_SOURCE_DIR}/brick.step")

add_executable(verisolid
	"main.cpp"
)

target_link_libraries(verisolid PUBLIC parasolid Eigen3::Eigen) #TKernel TKXSBase TKSTEP TKBRep TKG3d TKMath TKMesh TKTopAlgo TKShHealing)
#target_include_directories(debug PUBLIC ${OpenCASCADE_INCLUDE_DIR})
#target_link_directories(debug PUBLIC ${OpenCASCADE_LIBRARY_DIR})
