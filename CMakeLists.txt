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

# Prefer static zlib — avoids GLIBC symbol-versioning issues in conda/spack
# environments on clusters. Falls back to dynamic if no static lib is found.
set(ZLIB_USE_STATIC_LIBS ON)
find_package(ZLIB QUIET)
if(NOT ZLIB_FOUND)
    set(ZLIB_USE_STATIC_LIBS OFF)
    find_package(ZLIB REQUIRED)
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
target_link_libraries(tuna PRIVATE xxhash ZLIB::ZLIB Threads::Threads)
