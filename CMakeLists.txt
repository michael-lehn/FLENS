cmake_minimum_required(VERSION 3.5...3.31)

project(FLENS C CXX)

set(INCLUDE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/include" CACHE PATH "Include files install path")

add_subdirectory(cxxblas)
add_subdirectory(cxxlapack)
add_subdirectory(cxxstd)
add_subdirectory(flens)
add_subdirectory(external)
