#include <stdio.h>
#include <string.h>
#include <time.h>

#include <stdlib.h>

#include "stb_image.h"
#include "stb_image_resize.h"

#define USE_ENCODE_FUNCTIONS
#define USE_UINT16_FUNCTIONS

#include "icer.h"
#include "color_util.h"

const char compressed_filename[] = "./compressed.bin";
const char filename[] = "./boatcolor.512.bmp";

void rgb888_packed_to_yuv(uint16_t *y_channel, uint16_t *u_channel, uint16_t *v_channel, uint8_t *img, size_t image_w, size_t image_h, size_t rowstride) {
    int32_t r, g, b;
    uint8_t *pixel;

    uint16_t *output_y, *output_u, *output_v;
    for (size_t row = 0;row < image_h;row++) {
        pixel = img + 3 * rowstride * row;
        output_y = y_channel + rowstride * row;
        output_u = u_channel + rowstride * row;
        output_v = v_channel + rowstride * row;
        for (size_t col = 0;col < image_w;col++) {
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

int main() {
    const size_t out_w = 512;
    const size_t out_h = 512;
    const int stages = 4;
    const enum icer_filter_types filt = ICER_FILTER_A;
    const int segments = 10;

    const int datastream_size = 100000;

    int src_w, src_h, n;
    uint8_t *data;

    uint8_t *resized = malloc(out_w*out_h*3);
    uint16_t *compress[3];

    for (int chan = 0;chan <= 2;chan++) {
        compress[chan] = malloc(out_w*out_h*2);
    }

    int res = 0;
    clock_t begin, end;

    icer_init();

    printf("test compression code\n");

    printf("loading image: \"%s\"\n", filename);
    data = stbi_load(filename, &src_w, &src_h, &n, 3);
    if (data == NULL) {
        printf("invalid image\nexiting...\n");
        return 0;
    }

    printf("loaded image\nwidth    : %5d\nheight   : %5d\nchannels : %5d\n", src_w, src_h, n);

    printf("resizing image to width: %4llu, height: %4llu\n", out_w, out_h);
    res = stbir_resize_uint8(data, src_w, src_h, 0,
                             resized, out_w, out_h, 0,
                             3);
    if (res == 0) {
        printf("resize failed\nexiting...\n");
        return 0;
    }
    printf("resize complete\n");

    printf("converting to yuv\n");
    rgb888_packed_to_yuv(compress[0], compress[1], compress[2], resized, out_w, out_h, out_w);

    uint8_t *datastream = malloc(datastream_size*2+500);
    icer_output_data_buf_typedef output;
    icer_init_output_struct(&output, datastream, datastream_size*2, datastream_size);

    begin = clock();
    icer_compress_image_yuv_uint16(compress[0], compress[1], compress[2], out_w, out_h, stages, filt, segments, &output);
    end = clock();

    printf("compressed size %llu, time taken: %lf\n", output.size_used, (float)(end-begin)/CLOCKS_PER_SEC);

    FILE *ptr1;

    ptr1 = fopen(compressed_filename,"wb");
    size_t written = fwrite(output.rearrange_start, sizeof(output.rearrange_start[0]), output.size_used, ptr1);
    printf("written: %llu\n", written);
    fflush(ptr1);
    fclose(ptr1);

    printf("output saved\n");

    free(resized);
    free(datastream);
    for (int chan = 0;chan <= 2;chan++) free(compress[chan]);
    stbi_image_free(data);

    return 0;
}
