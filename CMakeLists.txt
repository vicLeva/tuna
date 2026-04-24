cmake_minimum_required(VERSION 3.17)
project(tuna VERSION 1.3.0 LANGUAGES CXX C)

# Default to Release if no build type is specified.
# Only set when tuna is the top-level project — do not override the parent's choice.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    endif()
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif()

# ── Dependencies ───────────────────────────────────────────────────────────────

find_package(Threads REQUIRED)
include(FetchContent)

# ── kff-cpp-api for KFF output format ─────────────────────────────────────────
# We skip kff's own CMakeLists (which would pull in unneeded benchmarks/executables)
# and build the library directly from the two source files we need.
FetchContent_Declare(
    kff_src
    GIT_REPOSITORY https://github.com/Kmer-File-Format/kff-cpp-api.git
    GIT_TAG        main
    GIT_SHALLOW    TRUE
)
FetchContent_GetProperties(kff_src)
if(NOT kff_src_POPULATED)
    FetchContent_Populate(kff_src)
endif()

set(KFF_VERSION_MAJOR 1)
set(KFF_VERSION_MINOR 0)
configure_file(
    ${kff_src_SOURCE_DIR}/kff_io.hpp.in
    ${CMAKE_BINARY_DIR}/kff_io.hpp
    @ONLY
)
add_library(kff_io STATIC ${kff_src_SOURCE_DIR}/kff_io.cpp)
target_compile_options(kff_io PRIVATE -w)
target_include_directories(kff_io SYSTEM PUBLIC
    ${kff_src_SOURCE_DIR}
    ${CMAKE_BINARY_DIR}
)

# ── zlib: prefer zlib-ng for ~2x faster gz decompression ──────────────────────
# zlib-ng compat mode is a drop-in replacement: same zlib.h API, same symbols,
# no code changes required. Disable with -DTUNA_USE_ZLIB_NG=OFF to fall back
# to the system zlib (e.g. if building in an environment where fetching is not
# possible or a system zlib-ng is already installed).
option(TUNA_USE_ZLIB_NG "Build against zlib-ng (faster gz) instead of system zlib" ON)

if(TUNA_USE_ZLIB_NG)
    include(FetchContent)
    FetchContent_Declare(
        zlib-ng
        GIT_REPOSITORY https://github.com/zlib-ng/zlib-ng.git
        GIT_TAG        2.2.2
        GIT_SHALLOW    TRUE
    )
    set(ZLIB_COMPAT       ON  CACHE INTERNAL "zlib-ng: enable zlib-compat API")
    set(ZLIB_ENABLE_TESTS OFF CACHE INTERNAL "zlib-ng: skip tests")
    set(BUILD_SHARED_LIBS OFF CACHE INTERNAL "zlib-ng: static only")
    FetchContent_MakeAvailable(zlib-ng)
    set(TUNA_ZLIB_TARGET zlibstatic)
else()
    # Fall back to system zlib. Prefer static to avoid GLIBC symbol-versioning
    # issues in conda/spack environments on clusters.
    set(ZLIB_USE_STATIC_LIBS ON)
    find_package(ZLIB QUIET)
    if(NOT ZLIB_FOUND)
        set(ZLIB_USE_STATIC_LIBS OFF)
        find_package(ZLIB REQUIRED)
    endif()
    set(TUNA_ZLIB_TARGET ZLIB::ZLIB)
endif()

# L1 cache line size required by kache-hash internals.
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    execute_process(
        COMMAND getconf LEVEL1_DCACHE_LINESIZE
        COMMAND tr -d "\n"
        OUTPUT_VARIABLE L1_CACHE_LINE_SIZE
    )
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    execute_process(
        COMMAND sysctl hw.cachelinesize
        COMMAND awk "{print $2}"
        COMMAND tr -d "\n"
        OUTPUT_VARIABLE L1_CACHE_LINE_SIZE
    )
endif()
# ── Compile flags ──────────────────────────────────────────────────────────────
#
# TUNA_SIMD_FLAGS   — instruction-set flags required to compile tuna's SIMD
#                     headers (helicase, kache-hash).  Propagated to users via
#                     INTERFACE because the headers themselves use the intrinsics.
#
# TUNA_PRIVATE_FLAGS — warnings, optimisations, constexpr limits.  Applied only
#                      to tuna's own compiled targets; never leaked to the parent.

if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64|i686")
    set(TUNA_SIMD_FLAGS -msse3 -msse4.1 -mavx -mavx2)
else()
    set(TUNA_SIMD_FLAGS "")
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CONSTEXPR_FLAG -fconstexpr-ops-limit=99000100)
else()
    set(CONSTEXPR_FLAG -fconstexpr-steps=99000100)
endif()
set(TUNA_PRIVATE_FLAGS -Wall -Wextra -Wno-unused-function ${CONSTEXPR_FLAG})

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    list(APPEND TUNA_PRIVATE_FLAGS -funroll-loops -march=native -fomit-frame-pointer)
endif()

# ── Targets ────────────────────────────────────────────────────────────────────

add_library(xxhash STATIC ${CMAKE_CURRENT_SOURCE_DIR}/external/xxHash/xxhash.c)
target_include_directories(xxhash PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_options(xxhash PRIVATE -O3)
set_target_properties(xxhash PROPERTIES C_STANDARD 11 C_STANDARD_REQUIRED ON)

# ── tuna::tuna — embeddable library target ─────────────────────────────────────
# Header-only INTERFACE target: all tuna logic lives in templates.
# Link against this to embed tuna's k-mer counting in your own project.
#
# Usage (CMake):
#   add_subdirectory(tuna)
#   target_link_libraries(my_target PRIVATE tuna::tuna)
#
# Usage (C++):
#   #include <tuna/tuna.hpp>
#   auto kmers = tuna::count_to<31>({"genome.fa"});
add_library(tuna_lib INTERFACE)
add_library(tuna::tuna ALIAS tuna_lib)
target_include_directories(tuna_lib INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/include/helicase_cpp
)
# C++ standard and definitions are propagated so users get the right dialect
# and tuna's required macros, but compile flags (warnings, optimisations) are
# not — those are tuna's internal concern.
target_compile_features(tuna_lib INTERFACE cxx_std_20)
target_compile_options(tuna_lib INTERFACE ${TUNA_SIMD_FLAGS})
target_compile_definitions(tuna_lib INTERFACE
    L1_CACHE_LINE_SIZE=${L1_CACHE_LINE_SIZE}
    TUNA_VERSION="${PROJECT_VERSION}"
    $<$<BOOL:${FIXED_K}>:FIXED_K=${FIXED_K}>
    $<$<BOOL:${FIXED_M}>:FIXED_M=${FIXED_M}>
)
target_link_libraries(tuna_lib INTERFACE xxhash ${TUNA_ZLIB_TARGET} Threads::Threads kff_io)

# ── tuna — CLI executable (only when tuna is the top-level project) ────────────
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    add_executable(tuna ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/CLI.cpp)
    target_compile_options(tuna PRIVATE ${TUNA_PRIVATE_FLAGS})
    target_link_libraries(tuna PRIVATE tuna_lib)

    # ── tuna-prot — protein k-mer counter ──────────────────────────────────────
    add_executable(tuna_prot
        ${CMAKE_CURRENT_SOURCE_DIR}/src/prot/main_prot.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/prot/CLI_prot.cpp)
    target_compile_features(tuna_prot PRIVATE cxx_std_20)
    target_compile_options(tuna_prot PRIVATE ${TUNA_PRIVATE_FLAGS})
    target_include_directories(tuna_prot PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/src/prot
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/include/helicase_cpp
    )
    target_compile_definitions(tuna_prot PRIVATE
        TUNA_VERSION="${PROJECT_VERSION}"
    )
    target_link_libraries(tuna_prot PRIVATE xxhash ${TUNA_ZLIB_TARGET} Threads::Threads)
endif()
