cmake_minimum_required(VERSION 3.23)
project(icer_compression C)

set(CMAKE_C_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_C_FLAGS "-Wall -Wextra")
set(CMAKE_C_FLAGS_DEBUG "-g")
set(CMAKE_C_FLAGS_RELEASE "-O2")

include_directories(example/inc)
include_directories(lib_icer/inc)
add_subdirectory(lib_icer)

add_executable(compress example/src/example_encode.c example/src/stb_lib.c example/inc/color_util.h)
add_executable(decompress example/src/example_decode.c example/src/stb_lib.c)

add_executable(compress_color example/src/example_encode_color.c example/src/stb_lib.c)
add_executable(decompress_color example/src/example_decode_color.c example/src/stb_lib.c)

target_link_libraries(compress PRIVATE m)
target_link_libraries(compress PUBLIC icer)

target_link_libraries(decompress PRIVATE m)
target_link_libraries(decompress PUBLIC icer)

target_link_libraries(compress_color PRIVATE m)
target_link_libraries(compress_color PUBLIC icer)

target_link_libraries(decompress_color PRIVATE m)
target_link_libraries(decompress_color PUBLIC icer)