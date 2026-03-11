cmake_minimum_required(VERSION 3.17)
project(tuna VERSION 1.0.0 LANGUAGES CXX C)

# ── Language standards ─────────────────────────────────────────────────────────
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# ── Build type ─────────────────────────────────────────────────────────────────
# Always pass -DCMAKE_BUILD_TYPE=Release explicitly. A stale CMakeCache.txt from
# a previous Debug configuration will otherwise produce an unoptimized binary
# (~6× slower). Plain `cmake ..` inherits the cached type.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ── Required system libraries ──────────────────────────────────────────────────
find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)

# ── L1 cache line size ────────────────────────────────────────────────────────
# Required by utility.hpp (Padded<T>) and kache-hash internals.
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

# ── Optional preprocessor defines ─────────────────────────────────────────────
# Number of k values instantiated at compile time (default: 32 instances).
if(INSTANCE_COUNT)
    add_compile_definitions(INSTANCE_COUNT=${INSTANCE_COUNT})
endif()

if(FIXED_K)
    add_compile_definitions(FIXED_K=${FIXED_K})
endif()

# ── Compile flags ──────────────────────────────────────────────────────────────
set(TUNA_COMPILE_FLAGS
    -Wall -Wextra
    -Wno-unused-function
    -funroll-loops -march=native -mavx -msse4.1 -msse3
    -fconstexpr-ops-limit=99000100
)

# ── xxHash (compiled C library; required by kache-hash) ───────────────────────
add_library(xxhash STATIC external/xxHash/xxhash.c)
target_include_directories(xxhash PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_compile_options(xxhash PRIVATE -O3)

# ── tuna binary ───────────────────────────────────────────────────────────────
add_executable(tuna src/main.cpp src/CLI.cpp)

target_include_directories(tuna PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/external/helicase_cpp
)

target_compile_options(tuna PRIVATE ${TUNA_COMPILE_FLAGS})

target_link_libraries(tuna PRIVATE
    xxhash
    ZLIB::ZLIB
    Threads::Threads
)
