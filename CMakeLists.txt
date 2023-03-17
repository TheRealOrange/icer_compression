cmake_minimum_required(VERSION 3.23)
project(icer_compression C)

set(CMAKE_C_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

include_directories( inc )

add_executable(icer_compression src/main.c inc/icer.h src/icer.c inc/stb_image_write.h inc/stb_image.h src/stb_lib.c inc/stb_image_resize.h)
