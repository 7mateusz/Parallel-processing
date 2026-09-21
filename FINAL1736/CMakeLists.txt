cmake_minimum_required(VERSION 3.16)
project(FINAL1736 CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_compile_options(-O3 -march=native -Wall -Wextra)

add_executable(sito-seq src/seq.cpp)
add_executable(sito-omp src/omp.cpp)
find_package(OpenMP REQUIRED)
target_link_libraries(sito-omp PRIVATE OpenMP::OpenMP_CXX)
