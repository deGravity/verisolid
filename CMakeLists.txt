cmake_minimum_required(VERSION 3.8)
project(verisolid)

set (CMAKE_CXX_STANDARD 11)

include(FetchContent)
FetchContent_Declare(
  ps_frustrum
  GIT_REPOSITORY   https://github.com/deGravity/parasolid_frustrum.git
  GIT_TAG   v1.0
)

FetchContent_MakeAvailable(ps_frustrum)


add_definitions(-DTEST_X_T="${CMAKE_CURRENT_SOURCE_DIR}/brick.x_t")

add_executable(verisolid "main.cpp")

target_link_libraries(verisolid PUBLIC parasolid_frustrum)
