#include <stdio.h>
#include <string.h>
#include <time.h>

#include <stdlib.h>

#include "stb_image_resize.h"
#include "stb_image_write.h"

#define USE_DECODE_FUNCTIONS
#define USE_UINT16_FUNCTIONS

#include "icer.h"

const char compressed_filename[] = "./compressed.bin";
const char filename[] = "./decompress.bmp";

int main() {
    const size_t out_w = 1000;
    const size_t out_h = 1000;
    const int stages = 4;
    const enum icer_filter_types filt = ICER_FILTER_A;
    const int segments = 6;

    uint16_t *decompress = malloc(out_w*out_h*2);
    uint8_t *display = malloc(out_w*out_h);

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
    res = icer_decompress_image_uint16(decompress, &decomp_w, &decomp_h, out_w*out_h, buf, length, stages, filt, segments);
    if (res != ICER_RESULT_OK) {
        printf("error: %d\n", res);
        return 0;
    }
    end = clock();
    printf("decompress time taken: %lf\n", (float)(end-begin)/CLOCKS_PER_SEC);

    for (size_t i = 0;i < out_w*out_h;i++) {
        display[i] = icer_min_int(0xff, decompress[i]);
    }

    printf("saving decompressed image to: \"%s\"\n", filename);
    res = stbi_write_bmp(filename, decomp_w, decomp_h, 1, display);
    if (res == 0) {
        printf("save failed\nexiting...\n");
        return 0;
    }

    free(decompress);
    free(display);
    free(buf);

    return 0;
}
