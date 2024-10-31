# ICER Image Compression Algorithm

The code in this repository implements the NASA ICER image compression algorithm as a C library. Said compression algorithm is a progressive, wavelet-based image compression algorithm designed to be resistant to data loss, making it suitable for use as the image compression algorithm when encoding images to be transmitted over unreliable delivery channels, such as those in satellite radio communications.

This library was designed with memory-constrained embedded systems in mind, hence the language choice of C, but it should function just as well on normal machines.

## Features

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

## How do I use it?

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

## How do I test it?

To run the example code, simply build it and place the `boat.512.bmp`  (and `boatcolor.512.bmp` for color encode) image as the same folder as the executable to generate the output.



### Note: Work in progress!

The library now supports:

- 8-bit and 16-bit image processing
- Full color support using YUV color space
- Multiple wavelet filter types
- Configurable error containment segments
- Progressive compression with size targeting

Future improvements may include:

- Additional color space support
- Enhanced compression options for different color channels

