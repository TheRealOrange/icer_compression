#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <stdint.h>
#include <getopt.h>

#include "stb_image.h"
#include "stb_image_write.h"
#include "color_util.h"

// ICER library configuration
#define USE_ENCODE_FUNCTIONS
#define USE_DECODE_FUNCTIONS
#define USE_UINT16_FUNCTIONS

#include "icer.h"

typedef struct {
  char *input_file;
  char *output_file;
  char *operation; // "compress" or "decompress"
  int stages;
  enum icer_filter_types filter;
  int segments;
  int target_size; // Target output size in bytes (0 = lossless)
  int force_color;
  int force_grayscale;
} icer_config_t;

void print_usage(const char *program_name) {
  printf("usage: %s <operation> <input> <output> [options]\n\n", program_name);
  printf("operations:\n");
  printf("  compress    Compress an image file\n");
  printf("  decompress  Decompress a compressed file\n\n");
  printf("compression options:\n");
  printf("  -s, --stages <n>      Number of wavelet decomposition stages (default: 4)\n");
  printf("  -f, --filter <type>   Filter type: A, B, C, D, E, F, Q (default: A)\n");
  printf("  -g, --segments <n>    Number of error containment segments (default: 6)\n");
  printf("  -t, --size <size>     Target compressed size in bytes (default: lossless)\n");
  printf("  -c, --color           Force color compression (YUV)\n");
  printf("  -G, --grayscale       Force grayscale compression\n");
  printf("\n");
  printf("decompression options:\n");
  printf("  -c, --color           Decompress as color image\n");
  printf("  -G, --grayscale       Decompress as grayscale image\n");
  printf("\n");
  printf("others:\n");
  printf("  --help                Show this help message\n\n");
  printf("For decompression, you must specify either --color or --grayscale\n");
}

enum icer_filter_types parse_filter_type(const char *filter_str) {
  if (strcasecmp(filter_str, "A") == 0) return ICER_FILTER_A;
  if (strcasecmp(filter_str, "B") == 0) return ICER_FILTER_B;
  if (strcasecmp(filter_str, "C") == 0) return ICER_FILTER_C;
  if (strcasecmp(filter_str, "D") == 0) return ICER_FILTER_D;
  if (strcasecmp(filter_str, "E") == 0) return ICER_FILTER_E;
  if (strcasecmp(filter_str, "F") == 0) return ICER_FILTER_F;
  if (strcasecmp(filter_str, "Q") == 0) return ICER_FILTER_Q;

  fprintf(stderr, "Invalid filter type: %s. Using default filter A.\n", filter_str);
  return ICER_FILTER_A;
}

void rgb888_packed_to_yuv(uint16_t *y_channel, uint16_t *u_channel, uint16_t *v_channel,
                          uint8_t *img, size_t image_w, size_t image_h, size_t rowstride) {
  int32_t r, g, b;
  uint8_t *pixel;
  uint16_t *output_y, *output_u, *output_v;

  for (size_t row = 0; row < image_h; row++) {
    pixel = img + 3 * rowstride * row;
    output_y = y_channel + rowstride * row;
    output_u = u_channel + rowstride * row;
    output_v = v_channel + rowstride * row;

    for (size_t col = 0; col < image_w; col++) {
      r = pixel[0];
      g = pixel[1];
      b = pixel[2];

      *output_y = CRGB2Y(r, g, b);
      *output_u = CRGB2Cb(r, g, b);
      *output_v = CRGB2Cr(r, g, b);

      pixel += 3;
      output_y++; output_u++; output_v++;
    }
  }
}

void yuv_to_rgb888_packed(uint16_t *y_channel, uint16_t *u_channel, uint16_t *v_channel,
                          uint8_t *img, size_t image_w, size_t image_h, size_t rowstride) {
  int32_t y, u, v;
  uint8_t *pixel;
  uint16_t *input_y, *input_u, *input_v;

  for (size_t row = 0; row < image_h; row++) {
    pixel = img + 3 * rowstride * row;
    input_y = y_channel + rowstride * row;
    input_u = u_channel + rowstride * row;
    input_v = v_channel + rowstride * row;

    for (size_t col = 0; col < image_w; col++) {
      y = *input_y;
      u = *input_u;
      v = *input_v;

      pixel[0] = CYCbCr2R(y, u, v);
      pixel[1] = CYCbCr2G(y, u, v);
      pixel[2] = CYCbCr2B(y, u, v);

      pixel += 3;
      input_y++; input_u++; input_v++;
    }
  }
}

int compress_image(icer_config_t *config) {
  int src_w, src_h, channels;
  uint8_t *data;
  clock_t begin, end;
  int result = 0;

  int target_channels = 0;
  if (config->force_color) {
    target_channels = 3;
  } else if (config->force_grayscale) {
    target_channels = 1;
  }

  data = stbi_load(config->input_file, &src_w, &src_h, &channels, target_channels);
  if (data == NULL) {
    fprintf(stderr, "Error: Could not load image %s\n", config->input_file);
    return 1;
  }

  printf("Loaded image: %s (%dx%d, %d channels)\n", config->input_file,  src_w, src_h, channels);

  int use_color = 0;
  if (target_channels == 0 && channels == 3) {
    use_color = 1;
  } else if (target_channels == 3) {
    use_color = 1;
  }

  printf("Compression mode: %s\n", use_color ? "Color (YUV)" : "Grayscale");

  uint16_t *compress[3] = {NULL, NULL, NULL};

  if (use_color) {
    for (int i = 0; i < 3; i++) {
      compress[i] = malloc(src_w * src_h * sizeof(uint16_t));
    }

    // Convert RGB to YUV for ICER
    rgb888_packed_to_yuv(compress[0], compress[1], compress[2], data, src_w, src_h, src_w);
  } else {
    compress[0] = malloc(src_w * src_h * sizeof(uint16_t));

    // convert to uint16 for ICER
    for (int i = 0; i < src_w * src_h; i++) {
      compress[0][i] = data[i];
    }
  }

  // calculate target size
  int byte_quote;
  if (config->target_size > 0) {
    // use specified target size
    byte_quote = config->target_size;
  } else {
    // lossless, set the quota the original image size to be safe
    byte_quote = src_w * src_h * (use_color ? 3 : 1);
  }

  // alloc output buffer to at least twice the target size
  int buffer_size = byte_quote * 2 + 50;
  uint8_t *datastream = malloc(buffer_size);

  icer_output_data_buf_typedef output;
  icer_init_output_struct(&output, datastream, buffer_size, byte_quote);

  printf("Starting compression...\n");
  printf("Parameters: stages=%d, filter=%d, segments=%d",
         config->stages, config->filter, config->segments);
  if (config->target_size > 0) {
    printf(", target_size=%.2fKB", config->target_size / 1024.0);
  } else {
    printf(", mode=lossless, quota=%.2fKB", byte_quote / 1024.0);
  }
  printf("\n");

  begin = clock();
  int compress_result;
  if (use_color) {
    compress_result = icer_compress_image_yuv_uint16(compress[0], compress[1], compress[2],
                                                     src_w, src_h, config->stages,
                                                     config->filter, config->segments, &output);
  } else {
    compress_result = icer_compress_image_uint16(compress[0], src_w, src_h, config->stages,
                                                 config->filter, config->segments, &output);
  }
  end = clock();

  if (compress_result != ICER_RESULT_OK && compress_result != ICER_BYTE_QUOTA_EXCEEDED) {
    fprintf(stderr, "Error: Compression failed with code %d\n", compress_result);
    result = 1;
    goto cleanup;
  }

  printf("Compression completed in %.3f seconds\n", (double)(end-begin)/CLOCKS_PER_SEC);
  printf("Compressed size: %zu bytes (%.1f%% of original)\n",
         output.size_used, 100.0 * output.size_used / (src_w * src_h * (use_color ? 3 : 1)));

  FILE *output_file = fopen(config->output_file, "wb");
  if (!output_file) {
    fprintf(stderr, "Error: Could not create output file %s\n", config->output_file);
    result = 1;
    goto cleanup;
  }

  size_t written = fwrite(output.rearrange_start, 1, output.size_used, output_file);
  fclose(output_file);

  if (written != output.size_used) {
    fprintf(stderr, "Error: Failed to write complete output file\n");
    result = 1;
    goto cleanup;
  }

  printf("Compressed image saved to: %s (%zu bytes)\n", config->output_file, output.size_used);

  cleanup:
  if (data) stbi_image_free(data);
  if (datastream) free(datastream);
  for (int i = 0; i < 3; i++) {
    if (compress[i]) free(compress[i]);
  }

  return result;
}

int decompress_image(icer_config_t *config) {
  FILE *input_file;
  uint8_t *compressed_data = NULL;
  size_t compressed_size = 0;
  size_t buffer_size = 1024;
  clock_t begin, end;
  int result = 0;

  if (!config->force_color && !config->force_grayscale) {
    fprintf(stderr, "Error: For decompression, you must specify either --color or --grayscale\n");
    return 1;
  }

  int use_color = config->force_color;
  input_file = fopen(config->input_file, "rb");
  if (!input_file) {
    fprintf(stderr, "Error: Could not open input file %s\n", config->input_file);
    return 1;
  }

  compressed_data = malloc(buffer_size);
  while (fread(compressed_data + compressed_size, 1, 1, input_file) == 1) {
    compressed_size++;
    if (compressed_size >= buffer_size - 1) {
      buffer_size += 1024;
      compressed_data = realloc(compressed_data, buffer_size);
    }
  }
  fclose(input_file);

  printf("Loaded compressed data: %zu bytes\n", compressed_size);

  uint16_t *decompress[3];
  uint8_t *display = NULL;

  // Get image dimensions
  size_t image_w, image_h;
  int dim_result = icer_get_image_dimensions(compressed_data, compressed_size, &image_w, &image_h);
  if (dim_result != ICER_RESULT_OK) {
    fprintf(stderr, "Error: Could not read image dimensions from compressed data\n");
    result = 1;
    goto cleanup;
  }

  printf("Image dimensions: %zux%zu\n", image_w, image_h);

  // Allocate decompression buffers
  for (int i = 0; i < 3; i++) {
    decompress[i] = malloc(image_w * image_h * sizeof(uint16_t));
  }

  printf("Starting decompression as %s image...\n", use_color ? "color" : "grayscale");

  begin = clock();

  size_t actual_w, actual_h;
  int decompress_result;

  if (use_color) {
    // Decompress as color
    decompress_result = icer_decompress_image_yuv_uint16(decompress[0], decompress[1], decompress[2],
                                                         &actual_w, &actual_h, image_w * image_h,
                                                         compressed_data, compressed_size);
    if (decompress_result == ICER_RESULT_OK) {
      display = malloc(image_w * image_h * 3);
      yuv_to_rgb888_packed(decompress[0], decompress[1], decompress[2], display, actual_w, actual_h, actual_w);
    }
  } else {
    // Decompress as grayscale
    decompress_result = icer_decompress_image_uint16(decompress[0], &actual_w, &actual_h, image_w * image_h,
                                                     compressed_data, compressed_size);
    if (decompress_result == ICER_RESULT_OK) {
      display = malloc(image_w * image_h);
      for (size_t i = 0; i < actual_w * actual_h; i++) {
        display[i] = (decompress[0][i] > 255) ? 255 : decompress[0][i];
      }
    }
  }

  end = clock();

  if (decompress_result != ICER_RESULT_OK) {
    fprintf(stderr, "Error: Decompression failed with code %d\n", decompress_result);
    fprintf(stderr, "This may indicate the compressed data was created in a different mode (%s vs %s)\n",
            use_color ? "color" : "grayscale", use_color ? "grayscale" : "color");
    result = 1;
    goto cleanup;
  }

  printf("Decompression completed in %.3f seconds\n", (double)(end-begin)/CLOCKS_PER_SEC);
  printf("Decompressed to %zux%zu %s image\n", actual_w, actual_h, use_color ? "color" : "grayscale");

  // Save decompressed image
  int save_result;
  if (use_color) {
    save_result = stbi_write_bmp(config->output_file, actual_w, actual_h, 3, display);
  } else {
    save_result = stbi_write_bmp(config->output_file, actual_w, actual_h, 1, display);
  }

  if (!save_result) {
    fprintf(stderr, "Error: Failed to save output image %s\n", config->output_file);
    result = 1;
    goto cleanup;
  }

  printf("Decompressed image saved to: %s\n", config->output_file);

  cleanup:
  if (compressed_data) free(compressed_data);
  if (display) free(display);
  for (int i = 0; i < 3; i++) {
    if (decompress[i]) free(decompress[i]);
  }

  return result;
}

int main(int argc, char *argv[]) {
  icer_config_t config = {
      .input_file = NULL,
      .output_file = NULL,
      .operation = NULL,
      .stages = 4,
      .filter = ICER_FILTER_A,
      .segments = 6,
      .target_size = 0, // 0 = lossless (no size limit)
      .force_color = 0,
      .force_grayscale = 0,
  };

  static struct option long_options[] = {
      {"stages", required_argument, 0, 's'},
      {"filter", required_argument, 0, 'f'},
      {"segments", required_argument, 0, 'g'},
      {"size", required_argument, 0, 't'},
      {"color", no_argument, 0, 'c'},
      {"grayscale", no_argument, 0, 'G'},
      {"help", no_argument, 0, 0},
      {0, 0, 0, 0}
  };

  // First, we need to determine the operation to know which parameters to allow
  if (argc < 2) {
    fprintf(stderr, "Error: Missing operation\n");
    print_usage(argv[0]);
    return 1;
  }

  char *operation = argv[1];
  int is_compress = (strcmp(operation, "compress") == 0);
  int is_decompress = (strcmp(operation, "decompress") == 0);

  if (!is_compress && !is_decompress) {
    fprintf(stderr, "Error: Operation must be 'compress' or 'decompress'\n");
    print_usage(argv[0]);
    return 1;
  }

  int option_index = 0;
  int c;

  while ((c = getopt_long(argc, argv, "s:f:g:t:cG", long_options, &option_index)) != -1) {
    switch (c) {
      case 's':
        if (!is_compress) {
          fprintf(stderr, "Error: --stages parameter only allowed for compression\n");
          return 1;
        }
        config.stages = atoi(optarg);
        if (config.stages < 1 || config.stages > 6) {
          fprintf(stderr, "Error: Stages must be between 1 and 6\n");
          return 1;
        }
        break;
      case 'f':
        if (!is_compress) {
          fprintf(stderr, "Error: --filter parameter only allowed for compression\n");
          return 1;
        }
        config.filter = parse_filter_type(optarg);
        break;
      case 'g':
        if (!is_compress) {
          fprintf(stderr, "Error: --segments parameter only allowed for compression\n");
          return 1;
        }
        config.segments = atoi(optarg);
        if (config.segments < 1 || config.segments > 32) {
          fprintf(stderr, "Error: Segments must be between 1 and 32\n");
          return 1;
        }
        break;
      case 't':
        if (!is_compress) {
          fprintf(stderr, "Error: --size parameter only allowed for compression\n");
          return 1;
        }
        config.target_size = atoi(optarg);
        if (config.target_size < 0) {
          fprintf(stderr, "Error: Target size must be non-negative (0 = lossless)\n");
          return 1;
        }
        break;
      case 'c':
        config.force_color = 1;
        break;
      case 'G':
        config.force_grayscale = 1;
        break;
      case 0:
        if (strcmp(long_options[option_index].name, "help") == 0) {
          print_usage(argv[0]);
          return 0;
        }
        break;
      case '?':
        return 1;
      default:
        print_usage(argv[0]);
        return 1;
    }
  }

  if (config.force_color && config.force_grayscale) {
    fprintf(stderr, "Error: Cannot specify both --color and --grayscale\n");
    return 1;
  }

  if (optind + 2 >= argc) {
    fprintf(stderr, "Error: Missing required arguments\n");
    print_usage(argv[0]);
    return 1;
  }

  config.operation = argv[optind];
  config.input_file = argv[optind + 1];
  config.output_file = argv[optind + 2];

  // init ICER library
  int init_result = icer_init();
  if (init_result != ICER_RESULT_OK) {
    fprintf(stderr, "Error: Failed to initialize ICER library\n");
    return 1;
  }

  if (strcmp(config.operation, "compress") == 0) {
    return compress_image(&config);
  } else {
    return decompress_image(&config);
  }
}