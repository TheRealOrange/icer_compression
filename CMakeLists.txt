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

add_executable(icer_compression example/src/main.c example/inc/stb_image_write.h example/inc/stb_image.h example/src/stb_lib.c example/inc/stb_image_resize.h)

target_link_libraries(icer_compression PRIVATE m)
target_link_libraries(icer_compression PUBLIC icer)