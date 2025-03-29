cmake_minimum_required(VERSION 3.0)

project(OPENMESHER2D VERSION 1.0 LANGUAGES CXX)

if (WIN32)
    add_definitions(
        -DNOMINMAX
        -D_USE_MATH_DEFINES)
endif ()

find_package(OpenMesh REQUIRED)

add_subdirectory(3rdparty/predicates)
add_subdirectory(src/geometry)
add_subdirectory(src/mesh)
add_subdirectory(src/mesher)
add_subdirectory(test/mesher)
