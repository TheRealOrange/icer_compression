cmake_minimum_required(VERSION 3.23)
project(icer C)

set(CMAKE_C_STANDARD 11)

include_directories(inc)

add_library(icer inc/icer.h src/icer_wavelet.c inc/crc.h src/crc32.c src/icer_encoding.c src/icer_decoding.c src/icer_partition.c src/icer_util.c src/icer_context_modeller.c src/icer_init.c src/icer_config.c src/icer_compress.c)