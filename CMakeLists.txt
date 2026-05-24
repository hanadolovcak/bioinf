cmake_minimum_required(VERSION 3.16)
project(bamboo_filter CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_compile_options(-Wall -Wextra -Wpedantic)

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -fsanitize=address -fno-omit-frame-pointer")
set(CMAKE_EXE_LINKER_FLAGS_DEBUG "-fsanitize=address")

add_library(bamboo STATIC
    src/hash.cpp
    src/segment.cpp
    src/bamboo_filter.cpp
)
target_include_directories(bamboo PUBLIC src/)

add_executable(benchmark
    benchmark/benchmark_main.cpp
    benchmark/synthetic_gen.cpp
    benchmark/fasta_parser.cpp
    benchmark/kmer_extractor.cpp
    benchmark/memory_util.cpp
)
target_link_libraries(benchmark PRIVATE bamboo)
target_include_directories(benchmark PRIVATE benchmark/ src/)

enable_testing()

foreach(t segment bamboo kmer fasta)
    add_executable(test_${t} test/test_${t}.cpp)
    target_link_libraries(test_${t} PRIVATE bamboo)
    if(t STREQUAL "kmer" OR t STREQUAL "fasta")
        target_sources(test_${t} PRIVATE
            benchmark/kmer_extractor.cpp
            benchmark/fasta_parser.cpp
            benchmark/synthetic_gen.cpp
        )
        target_include_directories(test_${t} PRIVATE benchmark/)
    endif()
    add_test(NAME ${t} COMMAND test_${t})
endforeach()

# Optional adapter for the reference impl (wanghanchengchn/bamboofilters).
# Builds only when the repo is alongside this project and target is x86_64.
set(REFERENCE_DIR "${CMAKE_SOURCE_DIR}/../bamboofilters" CACHE PATH
    "Path to a checkout of wanghanchengchn/bamboofilters")

if(EXISTS "${REFERENCE_DIR}/src/bamboofilter/bamboofilter.hpp")
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
        message(STATUS "Reference adapter: ENABLED (REFERENCE_DIR=${REFERENCE_DIR})")
        add_executable(reference_adapter
            extern/reference_adapter.cpp
            ${REFERENCE_DIR}/src/common/BOBHash.cpp
            benchmark/synthetic_gen.cpp
            benchmark/fasta_parser.cpp
            benchmark/memory_util.cpp
        )
        target_include_directories(reference_adapter PRIVATE
            ${REFERENCE_DIR}/src
            benchmark/
        )
        target_compile_options(reference_adapter PRIVATE
            -mavx2 -fpermissive
            -Wno-unused-result -Wno-narrowing -Wno-pedantic
        )
    else()
        message(STATUS "Reference adapter: DISABLED (need x86_64, have ${CMAKE_SYSTEM_PROCESSOR})")
    endif()
else()
    message(STATUS "Reference adapter: DISABLED (${REFERENCE_DIR} not found)")
endif()
