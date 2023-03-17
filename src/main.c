#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>

#include "stb_image.h"
#include "stb_image_resize.h"
#include "stb_image_write.h"

#include "icer.h"

const char filename[] = "../lena_color.gif";
const char resized_filename[] = "../pic_resized.bmp";
const char resized_6bpp_filename[] = "../pic_resized_6bpp.bmp";
const char wavelet_filename[] = "../wavelet_img.bmp";
const char wavelet_inv_filename[] = "../wavelet_inv_img.bmp";

void reduce_bit_depth(unsigned char *buf, size_t len, uint8_t bits) {
    unsigned char *it = buf;
    for (;it < (buf+len);it++) {
        (*it) >>= bits;
    }
}

void increase_bit_depth(unsigned char *buf, size_t len, uint8_t bits) {
    unsigned char *it = buf;
    for (;it < (buf+len);it++) {
        (*it) <<= bits;
    }
}

bool compare(const unsigned char *buf1, const unsigned char *buf2, size_t len, size_t w) {
    bool identical = true;
    size_t cnt = 0;
    for (size_t it = 0;it < len;it++) {
        if (buf1[it] != buf2[it]) {
            size_t row = it / w;
            size_t col = it % w;
            //printf("diff at %4zu: src=%3d dst=%3d row=%3d col=%3d\n", it, buf1[it], buf2[it], row, col);
            cnt++;
            identical = false;
        }
    }
    printf("total differences: %zu\n", cnt);

    return identical;
}

int main() {
    const size_t out_w = 1000;
    const size_t out_h = 1500;
    const size_t out_channels = 1;

    int src_w, src_h, n;
    uint8_t *data;

    uint8_t *resized = malloc(out_w*out_h*out_channels);
    uint8_t *transformed = malloc(out_w*out_h*out_channels);

    int res = 0;
    clock_t begin, end;

    printf("test compression code\n");

    printf("loading image: \"%s\"\n", filename);
    data = stbi_load(filename, &src_w, &src_h, &n, out_channels);
    if (data == NULL) {
        printf("invalid image\nexiting...\n");
        return 0;
    }

    printf("loaded image\nwidth    : %5d\nheight   : %5d\nchannels : %5d\nout_chn  : %5d\n", src_w, src_h, n, out_channels);

    printf("resizing image to width: %4d, height: %4d\n", out_w, out_h);
    res = stbir_resize_uint8(data, src_w, src_h, 0,
                             resized, out_w, out_h, 0,
                             out_channels);
    if (res == 0) {
        printf("resize failed\nexiting...\n");
        return 0;
    }

    printf("saving resized image to: \"%s\"\n", resized_filename);
    res = stbi_write_bmp(resized_filename, out_w, out_h, out_channels, resized);
    if (res == 0) {
        printf("save failed\nexiting...\n");
        return 0;
    }
    int bit_red = 2;
    reduce_bit_depth(resized, out_w*out_h*out_channels, bit_red);
    increase_bit_depth(resized, out_w*out_h*out_channels, bit_red);

    printf("saving bit-depth reduced image to: \"%s\"\n", resized_6bpp_filename);
    res = stbi_write_bmp(resized_6bpp_filename, out_w, out_h, out_channels, resized);
    if (res == 0) {
        printf("save failed\nexiting...\n");
        return 0;
    }

    memcpy(transformed, resized, out_h*out_w*out_channels);

    reduce_bit_depth(transformed, out_w*out_h*out_channels, bit_red);

    bool overflow;
    begin = clock();
    overflow = icer_wavelet_transform_stages(transformed, out_w, out_h, 8, ICER_FILTER_F);
    end = clock();
    printf("time taken for wavelet transform: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);
    if (overflow) {
        printf("overflow\n");
    }

    uint8_t *datastart = transformed + icer_ceil_div_size_t(out_w, 2);
    partition_param_typdef partition;
    icer_generate_partition_parameters(&partition, icer_floor_div_size_t(out_w, 2), icer_ceil_div_size_t(out_h, 2), 20);
    for (int i = 0;i < 7;i++) {
        printf("lsb: %2d\n", i);
        compress_partition_uint8(datastart, &partition, out_w, ICER_SUBBAND_HL, i);
    }

    printf("saving wavelet transformed image to: \"%s\"\n", wavelet_filename);
    res = stbi_write_bmp(wavelet_filename, out_w, out_h, out_channels, transformed);
    if (res == 0) {
        printf("save failed\nexiting...\n");
        return 0;
    }

    begin = clock();
    overflow = icer_inverse_wavelet_transform_stages(transformed, out_w, out_h, 8, ICER_FILTER_F);
    end = clock();
    printf("time taken for inverse wavelet transform: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);
    if (overflow) {
        printf("overflow\n");
    }
    increase_bit_depth(transformed, out_w*out_h*out_channels, bit_red);

    printf("saving inverse wavelet transformed image to: \"%s\"\n", wavelet_inv_filename);
    res = stbi_write_bmp(wavelet_inv_filename, out_w, out_h, out_channels, transformed);
    if (res == 0) {
        printf("save failed\nexiting...\n");
        return 0;
    }


    if (compare(transformed, resized, out_w*out_h*out_channels, out_w)) {
        printf("result is identical\n");
    } else {
        printf("result is different\n");
    }

    free(resized);
    free(transformed);
    stbi_image_free(data);

    /*
    char str[1000], stri[1000];
    for (int i = 3;i < 60;i++) {
        int j = 0;
        for (;j < i;j++) {
            if (j < (i/2)+(i%2)) {
                str[j] = (char)((int)'a' + (j%26));
            } else {
                str[j] = (char)((int)'A' + ((j-((i/2)+(i%2)))%26));
            }

        }
        str[j] = '\0';
        printf("%d\n%s\n", i, str);
        strcpy(stri, str);
        icer_interleave_uint8((uint8_t*)stri, strlen(stri), 1);
        printf("%s\n", stri);
        icer_deinterleave_uint8((uint8_t*)stri, strlen(stri), 1);
        if (strcmp(str, stri) != 0) {
            printf("%s\n\n", stri);
        }
    }*/

    return 0;
}
