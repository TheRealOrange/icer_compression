#include <stdio.h>
#include <string.h>
#include <time.h>

#include <stdlib.h>

#include "stb_image_resize.h"
#include "stb_image_write.h"

#define USE_DECODE_FUNCTIONS
#define USE_UINT16_FUNCTIONS

#include "icer.h"
#include "color_util.h"

const char compressed_filename[] = "./compressed.bin";
const char filename[] = "./decompress.bmp";

void yuv_to_rgb888_packed(uint16_t *y_channel, uint16_t *u_channel, uint16_t *v_channel, uint8_t *img, size_t image_w, size_t image_h, size_t rowstride) {
    int32_t y, u, v;
    uint8_t *pixel;

    uint16_t *input_y, *input_u, *input_v;
    for (size_t row = 0;row < image_h;row++) {
        pixel = img + 3 * rowstride * row;
        input_y = y_channel + rowstride * row;
        input_u = u_channel + rowstride * row;
        input_v = v_channel + rowstride * row;
        for (size_t col = 0;col < image_w;col++) {
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

int main() {
    const size_t out_w = 1000;
    const size_t out_h = 1000;

    uint16_t *decompress[3];
    uint8_t *display = malloc(out_w*out_h*3);

    for (int chan = 0;chan <= 2;chan++) {
        decompress[chan] = malloc(out_w*out_h*2);
    }

    int res = 0;
    clock_t begin, end;

    icer_init();

    printf("test decompression code\n");

    FILE *ptr2;
    uint8_t *buf = malloc(500);
    size_t buf_size = 500;
    size_t length = 0;

    ptr2 = fopen(compressed_filename,"rb");
    while (fread(buf+length, sizeof *buf, 1, ptr2) == 1) {
        if (length >= buf_size-1) {
            buf_size += 500;
            buf = (uint8_t*)realloc(buf, buf_size);
        }
        length++;
    }
    fclose(ptr2);

    printf("decompress start\n");

    size_t decomp_w, decomp_h;
    begin = clock();
    res = icer_decompress_image_yuv_uint16(decompress[0], decompress[1], decompress[2], &decomp_w, &decomp_h, out_w*out_h, buf, length);
    end = clock();
    if (res != ICER_RESULT_OK) {
        printf("error: %d\n", res);
        return 0;
    }
    printf("decompress time taken: %lf\n", (float)(end-begin)/CLOCKS_PER_SEC);

    printf("converting to rgb\n");
    yuv_to_rgb888_packed(decompress[0], decompress[1], decompress[2], display, decomp_w, decomp_h, decomp_w);

    printf("saving decompressed image to: \"%s\"\n", filename);
    res = stbi_write_bmp(filename, decomp_w, decomp_h, 3, display);
    if (res == 0) {
        printf("save failed\nexiting...\n");
        return 0;
    }

    for (int chan = 0;chan <= 2;chan++) free(decompress[chan]);
    free(display);
    free(buf);

    return 0;
}
