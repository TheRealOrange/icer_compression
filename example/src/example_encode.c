#include <stdio.h>
#include <string.h>
#include <time.h>

#include <stdlib.h>

#include "stb_image.h"
#include "stb_image_resize.h"

#define USE_ENCODE_FUNCTIONS
#define USE_UINT16_FUNCTIONS

#include "icer.h"

const char compressed_filename[] = "./compressed.bin";
const char filename[] = "./boat.512.bmp";

int main() {
    const size_t out_w = 512;
    const size_t out_h = 512;
    const int stages = 4;
    const enum icer_filter_types filt = ICER_FILTER_A;
    const int segments = 6;

    const int datastream_size = 30000;

    int src_w, src_h, n;
    uint8_t *data;

    uint8_t *resized = malloc(out_w*out_h);
    uint16_t *compress = malloc(out_w*out_h*2);

    int res = 0;
    clock_t begin, end;

    icer_init();

    printf("test compression code\n");

    printf("loading image: \"%s\"\n", filename);
    data = stbi_load(filename, &src_w, &src_h, &n, 1);
    if (data == NULL) {
        printf("invalid image\nexiting...\n");
        return 0;
    }

    printf("loaded image\nwidth    : %5d\nheight   : %5d\nchannels : %5d\n", src_w, src_h, n);

    printf("resizing image to width: %4zu, height: %4zu\n", out_w, out_h);
    res = stbir_resize_uint8(data, src_w, src_h, 0,
                             resized, out_w, out_h, 0,
                             1);
    if (res == 0) {
        printf("resize failed\nexiting...\n");
        return 0;
    }
    printf("resize complete\n");

    printf("converting to int16\n");
    for (size_t i = 0;i < out_h*out_w;i++) {
        compress[i] = resized[i];
    }

    uint8_t *datastream = malloc(datastream_size*2+500);
    icer_output_data_buf_typedef output;
    icer_init_output_struct(&output, datastream, datastream_size*2, datastream_size);

    begin = clock();
    icer_compress_image_uint16(compress, out_w, out_h, stages, filt, segments, &output);
    end = clock();

    printf("compressed size %zu, time taken: %lf\n", output.size_used, (float)(end-begin)/CLOCKS_PER_SEC);

    FILE *ptr1;

    ptr1 = fopen(compressed_filename,"wb");
    size_t written = fwrite(output.rearrange_start, sizeof(output.rearrange_start[0]), output.size_used, ptr1);
    printf("written: %zu\n", written);
    fflush(ptr1);
    fclose(ptr1);

    printf("output saved\n");

    free(resized);
    free(compress);
    free(datastream);
    stbi_image_free(data);

    return 0;
}
