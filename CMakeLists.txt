cmake_minimum_required(VERSION 2.8.12)

project(SeismicMesh)

##############################################################################

# Add custom CMake modules
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# Some screen output
message(STATUS "OS detected: ${CMAKE_SYSTEM_NAME}")
message(STATUS "CXX Compiler detected: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMake additional search path for libraries: ${CMAKE_LIBRARY_PATH}")

# CGAL and its components
find_package(CGAL)
message(STATUS "CGAL version: ${CGAL_VERSION}")

if( CGAL_VERSION  VERSION_LESS 5.0)
    message(FATAL_ERROR "This project requires at least CGAL 5.0 library and will not be compiled.")
endif()


# include helper file
include( ${CGAL_USE_FILE})

set (MIGRATION_SRCE "SeismicMesh/migration/cpp")
set (GENERATION_SRCE "SeismicMesh/generation/cpp")
set (SIZING_SRCE "SeismicMesh/sizing/cpp")
set (GEOMETRY_SRCE "SeismicMesh/geometry/cpp")

include_directories (${MIGRATION_SRCE} ${GENERATION_SRCE} ${SIZING_SRC} ${GEOMETRY_SRCE})

find_package(pybind11)
message(STATUS "Found pybind11 v${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")

add_subdirectory(pybind11)
pybind11_add_module(cpputils ${SOURCES}  "${MIGRATION_SRCE}/cpputils.cpp")
pybind11_add_module(FastHJ ${SOURCES} "${SIZING_SRCE}/FastHJ.cpp")
pybind11_add_module(c_cgal ${SOURCES} "${GENERATION_SRCE}/delaunay.cpp")
pybind11_add_module(delaunay_class ${SOURCES} "${GENERATION_SRCE}/delaunay_class.cpp")
pybind11_add_module(delaunay_class3 ${SOURCES} "${GENERATION_SRCE}/delaunay_class3.cpp")
pybind11_add_module(fast_geometry ${SOURCES} "${GEOMETRY_SRCE}/fast_geometry.cpp")
