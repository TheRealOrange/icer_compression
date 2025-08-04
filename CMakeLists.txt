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

# Implementation examples
add_executable(compress example/src/example_encode.c example/src/stb_lib.c)
add_executable(decompress example/src/example_decode.c example/src/stb_lib.c)
add_executable(compress_color example/src/example_encode_color.c example/src/stb_lib.c)
add_executable(decompress_color example/src/example_decode_color.c example/src/stb_lib.c)

# Command-line utility
add_executable(icer_util example/src/icer_util.c example/src/stb_lib.c)

# Examples
target_link_libraries(compress PRIVATE m)
target_link_libraries(compress PUBLIC icer)

target_link_libraries(decompress PRIVATE m)
target_link_libraries(decompress PUBLIC icer)

target_link_libraries(compress_color PRIVATE m)
target_link_libraries(compress_color PUBLIC icer)

target_link_libraries(decompress_color PRIVATE m)
target_link_libraries(decompress_color PUBLIC icer)

# Utility
target_link_libraries(icer_util PRIVATE m)
target_link_libraries(icer_util PUBLIC icer)

#  testing with sample images
add_custom_target(test_examples
        COMMAND ${CMAKE_COMMAND} -E echo "Testing ICER implementation examples..."

        # Copy test images from source to build directory if they exist
        COMMAND test -f ${CMAKE_CURRENT_SOURCE_DIR}/boat.512.bmp && ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/boat.512.bmp . || echo "boat.512.bmp not found in source directory"
        COMMAND test -f ${CMAKE_CURRENT_SOURCE_DIR}/boatcolor.512.bmp && ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/boatcolor.512.bmp . || echo "boatcolor.512.bmp not found in source directory"

        COMMAND ${CMAKE_COMMAND} -E echo "Testing grayscale compression:"
        COMMAND test -f boat.512.bmp && ./compress && ./decompress || echo "boat.512.bmp not available, skipping grayscale test"
        COMMAND ${CMAKE_COMMAND} -E echo ""
        COMMAND ${CMAKE_COMMAND} -E echo "Testing color compression:"
        COMMAND test -f boatcolor.512.bmp && ./compress_color && ./decompress_color || echo "boatcolor.512.bmp not available, skipping color test"
        COMMAND ${CMAKE_COMMAND} -E echo ""
        COMMAND ${CMAKE_COMMAND} -E echo "Testing standalone utility:"
        COMMAND test -f boat.512.bmp && ./icer_util compress boat.512.bmp test_compressed.bin --grayscale && ./icer_util decompress test_compressed.bin test_decompressed.bmp --grayscale || echo "boat.512.bmp not available, skipping utility test"
        COMMAND test -f boatcolor.512.bmp && ./icer_util compress boatcolor.512.bmp test_compressed_color.bin --color && ./icer_util decompress test_compressed_color.bin test_decompressed_color.bmp --color || echo "boatcolor.512.bmp not available, skipping color utility test"
        DEPENDS compress decompress compress_color decompress_color icer_util
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Running ICER compression tests"
)