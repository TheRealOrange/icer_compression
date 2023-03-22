cmake_minimum_required(VERSION 3.23)
project(icer_compression C)

set(CMAKE_C_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_C_FLAGS "-Wall -Wextra")
set(CMAKE_C_FLAGS_DEBUG "-g")
set(CMAKE_C_FLAGS_RELEASE "-O2")

include_directories( inc )

add_executable(icer_compression src/main.c inc/icer.h src/icer_wavelet.c inc/stb_image_write.h inc/stb_image.h src/stb_lib.c inc/stb_image_resize.h inc/crc.h src/crc32.c src/icer_encoding.c src/icer_decoding.c src/icer_partition.c src/icer_util.c src/icer_context_modeller.c src/icer_init.c src/icer_config.c src/icer_compress.c)
