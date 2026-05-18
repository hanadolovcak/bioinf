cmake_minimum_required(VERSION 3.14)
project(bioinf)

set(CMAKE_CXX_STANDARD 17)

add_executable(bioinf
    src/main.cpp
    src/BambooFilter.cpp
    src/CuckooFilter.cpp
)

target_include_directories(bioinf PRIVATE src/)
