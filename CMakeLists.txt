cmake_minimum_required(VERSION 3.10)
set(CMAKE_GENERATOR "MinGW Makefiles")
message("using ${CMAKE_GENERATOR} for makefile generation")

# set the project name
project (LammpsToVTK)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp ${PROJECT_SOURCE_DIR}/src/*.hpp)

add_compile_options("$<$<CONFIG:RELEASE>:-W;-Wall;-Wfloat-conversion;-O3>")
add_compile_options("$<$<CONFIG:DEBUG>:-W;-Wall;-Wfloat-conversion;-O0;-g>")


add_executable(${PROJECT_NAME} ${sources})
