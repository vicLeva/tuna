cmake_minimum_required(VERSION 3.17)
project(tuna VERSION 1.0.0 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Default to Release if no build type is specified.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# ── Dependencies ───────────────────────────────────────────────────────────────

find_package(Threads REQUIRED)

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
add_compile_definitions(L1_CACHE_LINE_SIZE=${L1_CACHE_LINE_SIZE})

# ── Optional: single-k build (-DFIXED_K=31) ───────────────────────────────────
add_compile_definitions(TUNA_VERSION="${PROJECT_VERSION}")

if(FIXED_K)
    add_compile_definitions(FIXED_K=${FIXED_K})
endif()

# ── Compile flags ──────────────────────────────────────────────────────────────

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CONSTEXPR_FLAG -fconstexpr-ops-limit=99000100)
else()
    set(CONSTEXPR_FLAG -fconstexpr-steps=99000100)
endif()
set(TUNA_COMPILE_FLAGS -Wall -Wextra -Wno-unused-function ${CONSTEXPR_FLAG})

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    list(APPEND TUNA_COMPILE_FLAGS -funroll-loops -march=native -fomit-frame-pointer)
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64|i686")
        list(APPEND TUNA_COMPILE_FLAGS -mavx -msse4.1 -msse3)
    endif()
endif()

# ── Targets ────────────────────────────────────────────────────────────────────

add_library(xxhash STATIC external/xxHash/xxhash.c)
target_include_directories(xxhash PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_compile_options(xxhash PRIVATE -O3)

add_executable(tuna src/main.cpp src/CLI.cpp)
target_include_directories(tuna PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include/helicase_cpp
)
target_compile_options(tuna PRIVATE ${TUNA_COMPILE_FLAGS})
target_link_libraries(tuna PRIVATE xxhash ${TUNA_ZLIB_TARGET} Threads::Threads)
