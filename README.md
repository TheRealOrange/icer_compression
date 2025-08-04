# ICER Image Compression Algorithm

The code in this repository implements the NASA ICER image compression algorithm as a C library. Said compression algorithm is a progressive, wavelet-based image compression algorithm designed to be resistant to data loss, making it suitable for use as the image compression algorithm when encoding images to be transmitted over unreliable delivery channels, such as those in satellite radio communications.

This library was designed with memory-constrained embedded systems in mind, hence the language choice of C, but it should function just as well on normal machines.

## Features

Note: You may want to check out the alternate branch `include-metadata` for future use, where the compressed output includes the details of the compression parameters used.

### Core Capabilities

- Full color image compression using YUV color space
- Progressive image compression with configurable quality levels
- Wavelet-based compression using various filter types
- Error containment segments to limit data loss effects
- Configurable output size targeting
- Support for both lossless and lossy compression
- Integer-only arithmetic operations for embedded systems
- Memory-conscious implementation with no dynamic allocation

## Compression Examples

Below are sample images of the compression algorithm producing output images of varying quality depending on the space allocated for it. 10 error containment segments are used, and 9 bitplanes were compressed

### Color Compression

| ![original_color](https://github.com/TheRealOrange/icer_compression/blob/master/assets/original_color.bmp?raw=true) | ![140kb_quota_color](https://github.com/TheRealOrange/icer_compression/blob/master/assets/140kb_quota_color.bmp?raw=true) |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| original color image (512x512, 768 kilobyte)                 | 140 kilobyte quota                                           |
| ![100kb_quota_color](https://github.com/TheRealOrange/icer_compression/blob/master/assets/100kb_quota_color.bmp?raw=true) | ![70kb_quota_color](https://github.com/TheRealOrange/icer_compression/blob/master/assets/70kb_quota_color.bmp?raw=true) |
| 100 kilobyte quota                                           | 70 kilobyte quota                                            |

### Grayscale Compression

| ![original](https://github.com/TheRealOrange/icer_compression/blob/master/assets/original.bmp?raw=true) | ![70kb_quota](https://github.com/TheRealOrange/icer_compression/blob/master/assets/70kb_quota.bmp?raw=true) |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| original image (512x512, 262 kilobyte)                       | 70 kilobyte quota                                            |
| ![50kb_quota](https://github.com/TheRealOrange/icer_compression/blob/master/assets/50kb_quota.bmp?raw=true) | ![30kb_quota](https://github.com/TheRealOrange/icer_compression/blob/master/assets/30kb_quota.bmp?raw=true) |
| 50 kilobyte quota                                            | 30 kilobyte quota                                            |

The compression effectiveness depends greatly on the number of error containment segments chosen.

## What is this?

This compression algorithm is based on [a document published by NASA](https://ipnpr.jpl.nasa.gov/progress_report/42-155/155J.pdf). This repository implements the algorithm as described by the paper, with additional support for color images.

This algorithm is best described by the abstract of the document by NASA, which reads as follows:

> ICER is a progressive, wavelet-based image data compressor designed to meet the specialized needs of deep-space applications while achieving state-of-the-art compression effectiveness. ICER can provide lossless and lossy compression, and incorporates an error-containment scheme to limit the effects of data loss during transmission. The Mars Exploration Rover (MER) mission will rely primarily on a software implementation of ICER for image compression. This article describes ICER and the methods it uses to meet its goals, and explains the rationale behind the choice of methods. Performance results also are presented.

## Credits

This library uses code from:

- CRC32 checksum library from the [SNIPPETS C Source Code Archive](https://github.com/vonj/snippets.org)
- Image parsing libraries from [STB libraries](https://github.com/nothings/stb)

## Key Features

There are a few key design considerations which were taking into account when writing this library, which are as follows:

- Designed to only utilise integer arithmetic operations, as designed by NASA
- Avoids dynamic memory allocation in the library, so memory can be better managed in memory-constrained embedded systems
- Optimised for execution on embedded systems
- Error containment features prevent data loss from affecting the quality of the whole image
- Ability to set a output size target and stop compression once the quota is reached
- Able to compress images losslessly (if the output quota permits)

(again, recommended that you carefully read the NASA document to understand why and how these are significant)

## How do I test it?

The project uses CMake for building. To build all components:

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

This will build:
- `compress`, `decompress`,  Grayscale compression/decompression **implementation** examples
- `compress_color`, `decompress_color` - Color compression/decompression **implementation** examples
- `icer_util` - Standalone command-line utility to compress/decompress ICER images

### Command-line Tool

The `icer_util` program is a command-line program for compressing and decompressing images with the ICER algorithm. You can use this to test the compression efficacy on input images, and play around with byte quotas.

#### Basic Usage

```bash
# Compress an image (with default parameters)
./icer_util compress input.jpg output.bin

# Decompress an image (compulsory to specify color or grayscale)
./icer_util decompress compressed.bin output.bmp --color
./icer_util decompress compressed.bin output.bmp --grayscale

# Show availabe options
./icer_util --help
```

#### Advanced Options

```bash
# Compress with specific parameters
# if compressing with non-default parameters it is
# necessary to supply the parameters used during decompression
./icer_util compress input.png output.bin \
    --stages 5 \
    --filter B \
    --segments 10 \
    -t 150000 # byte quota
    
# Decompress with specific parameters
# if the image was compressed with non-default parameters it is
# necessary to supply the parameters used
./icer_util decompress output.bin output.bmp \
    --stages 5 \
    --filter B \
    --segments 10

# Set target compression size in bytes
./icer_util compress input.jpg output.bin -t 100000
```

#### Supported Image Formats

- **Input (compression):** JPEG, PNG, BMP, TGA, and other formats supported by the STB Image library
- **Output (decompression):** BMP format

### Integration Examples

The integration examples show an example usage of the library functions within a project. It is for reference only, though you can build and run it to test.

To run the example code, simply build it and place the `boat.512.bmp`  (and `boatcolor.512.bmp` for color encode) image as the same folder as the executable to generate the output.

The examples demonstrate:
- `compress` / `example_encode.c` - Basic grayscale image compression
- `decompress` / `example_decode.c` - Basic grayscale image decompression
- `compress_color` / `example_encode_color.c` - Color image compression using YUV
- `decompress_color` / `example_decode_color.c` - Color image decompression from YUV

## How do I use it in my project?

The project is designed as a C-library and as such is intended to be used as a part of a larger program, not as a standalone program. 

To include the library, select the appropriate flags for which parts of the library to compile:

```c
#define USE_ENCODE_FUNCTIONS // Compile encode functions
#define USE_DECODE_FUNCTIONS // Compile decode functions

#define USE_UINT16_FUNCTIONS // Compile functions for uint16_t arrays
#define USE_UINT8_FUNCTIONS  // Compile functions for uint8_t arrays

#include "icer.h"
```

Then, somewhere at the start, one must call 

```c
int icer_init();
```

in order to initialise the constants and look up tables used in encode and decode functions. 

The main functions which are important for the usage of the library are (here we only show the `uint8_t` functions, but the `uint16_t` versions exists too):

1. Compression functions

```c
// For grayscale/single channel images
int icer_compress_image_uint8(uint8_t *image, size_t image_w, size_t image_h, uint8_t stages, enum icer_filter_types filt, uint8_t segments, icer_output_data_buf_typedef *output_data);

// For colour images
int icer_compress_image_yuv_uint8(uint8_t *y_channel, uint8_t *u_channel, uint8_t *v_channel, size_t image_w, size_t image_h, uint8_t stages, enum icer_filter_types filt, uint8_t segments, icer_output_data_buf_typedef *output_data);
```

which is the function which enables the user to specify an image buffer of a specific width and height, as well as the filter coefficients to use, the number of time to perform wavelet decomposition, and the number of error containment segments to subdivide the image into

2. Decompression functions

```c
// For grayscale/single channel images
int icer_decompress_image_uint8(uint8_t *image, size_t *image_w, size_t *image_h, size_t image_bufsize, uint8_t *datastream, size_t data_length, uint8_t stages, enum icer_filter_types filt, uint8_t segments);

// For colour images
int icer_decompress_image_yuv_uint8(uint8_t *y_channel, uint8_t *u_channel, uint8_t *v_channel, size_t *image_w, size_t *image_h, size_t image_bufsize, const uint8_t *datastream, size_t data_length, uint8_t stages, enum icer_filter_types filt, uint8_t segments);
```

which is the function which enables the user to decompress a byte stream stored inside `datastream` and fill the image and the dimensions into the buffer specified.

3. Utility Functions

```c
int icer_get_image_dimensions(const uint8_t *datastream, size_t data_length, size_t *image_w, size_t *image_h);
```

which can be used to get the dimensions of an encoded image (useful for knowing how much space to allocate for the image to decode it).

## Advanced Configuration

The library provides several compile-time configuration options through preprocessor definitions that allow you to optimize memory usage and customize functionality for your specific needs.

### Memory and Buffer Configuration

```c
// Maximum number of error containment segments
#ifndef ICER_MAX_SEGMENTS
#define ICER_MAX_SEGMENTS 32
#endif

// Maximum number of wavelet decomposition stages
#ifndef ICER_MAX_DECOMP_STAGES
#define ICER_MAX_DECOMP_STAGES 6
#endif

// Maximum number of packets for 8-bit and 16-bit processing
#ifndef ICER_MAX_PACKETS
#define ICER_MAX_PACKETS 300
#endif

#ifndef ICER_MAX_PACKETS_16
#define ICER_MAX_PACKETS_16 800
#endif

// Number of bitplanes to compress for 8-bit and 16-bit modes
#ifndef ICER_BITPLANES_TO_COMPRESS_8
#define ICER_BITPLANES_TO_COMPRESS_8 7
#endif

#ifndef ICER_BITPLANES_TO_COMPRESS_16
#define ICER_BITPLANES_TO_COMPRESS_16 9
#endif
```

### Feature Selection Flags

```c
// Select which bit depth functions to compile
#if !defined(USE_UINT8_FUNCTIONS) && !defined(USE_UINT16_FUNCTIONS)
#define USE_UINT8_FUNCTIONS
#define USE_UINT16_FUNCTIONS
#endif

// Select whether to compile encode/decode functions
#if !defined(USE_DECODE_FUNCTIONS) && !defined(USE_ENCODE_FUNCTIONS)
#define USE_DECODE_FUNCTIONS
#define USE_ENCODE_FUNCTIONS
#endif
```

### Advanced Memory Management

```c
// Enable user-provided buffers for fine-grained memory control
//#define USER_PROVIDED_BUFFERS
```

When `USER_PROVIDED_BUFFERS` is defined, you must allocate the following buffers:

For 8-bit encoding:

```c
icer_packet_context icer_packets[ICER_MAX_PACKETS];
icer_image_segment_typedef *icer_rearrange_segments_8[ICER_CHANNEL_MAX + 1][ICER_MAX_DECOMP_STAGES + 1][ICER_SUBBAND_MAX + 1][7][ICER_MAX_SEGMENTS + 1];
```

For 8-bit decoding:

```c
icer_image_segment_typedef *icer_reconstruct_data_8[ICER_CHANNEL_MAX + 1][ICER_MAX_DECOMP_STAGES + 1][ICER_SUBBAND_MAX + 1][ICER_MAX_SEGMENTS + 1][7];
```

For 16-bit encoding:

```c
icer_packet_context icer_packets_16[ICER_MAX_PACKETS_16];
icer_image_segment_typedef *icer_rearrange_segments_16[ICER_CHANNEL_MAX + 1][ICER_MAX_DECOMP_STAGES + 1][ICER_SUBBAND_MAX + 1][15][ICER_MAX_SEGMENTS + 1];
```

For 16-bit decoding:

```c
icer_image_segment_typedef *icer_reconstruct_data_16[ICER_CHANNEL_MAX + 1][ICER_MAX_DECOMP_STAGES + 1][ICER_SUBBAND_MAX + 1][ICER_MAX_SEGMENTS + 1][15];
```

Common encoding buffer:

```c
uint16_t icer_encode_circ_buf[ICER_CIRC_BUF_SIZE];  // ICER_CIRC_BUF_SIZE = 2048
```

### Custom CRC32 Implementation

```c
// Define custom CRC32 implementation (e.g., hardware acceleration)
//#define CRC32BUF_FUNCTION(x, y) your_custom_crc32_function
```

When defined, you can provide your own CRC32 implementation with the following signature:

```c
uint32_t crc32buf(char *buf, size_t len)
```

## Note: Work in progress!

The library now supports:

- 8-bit and 16-bit image processing
- Full color support using YUV color space
- Multiple wavelet filter types
- Configurable error containment segments
- Progressive compression with size targeting

Future improvements may include:

- Additional color space support
- Enhanced compression options for different color channels

