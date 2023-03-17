#ifndef ICER_COMPRESSION_ICER_H
#define ICER_COMPRESSION_ICER_H

#include <stdbool.h>
#include <stdlib.h>

#define MAX_K 12

enum icer_status {
    ICER_RESULT_OK = 0,
    ICER_INTEGER_OVERFLOW = 1,
    ICER_SIZE_ERROR = 2,
    ICER_TOO_MANY_SEGMENTS = 3,
    ICER_TOO_MANY_STAGES = 4,
    ICER_CACHE_TOO_SMALL = 5
};

enum icer_filter_types {
    ICER_FILTER_A = 0,
    ICER_FILTER_B,
    ICER_FILTER_C,
    ICER_FILTER_D,
    ICER_FILTER_E,
    ICER_FILTER_F,
    ICER_FILTER_Q
};

enum icer_filter_params {
    ICER_FILTER_COEF_ALPHA_N1 = 0,
    ICER_FILTER_COEF_ALPHA_0,
    ICER_FILTER_COEF_ALPHA_1,
    ICER_FILTER_COEF_BETA
};

#define ICER_MAX_SEGMENTS 32
#define ICER_FILTER_DENOMINATOR 16

typedef struct {
    uint16_t w;
    uint16_t h;
    uint16_t r;
    uint16_t c;
    uint16_t r_t;
    uint16_t h_t;
    uint16_t x_t;
    uint16_t c_t0;
    uint16_t y_t;
    uint16_t r_t0;
    uint16_t x_b;
    uint16_t c_b0;
    uint16_t y_b;
    uint16_t r_b0;
    uint16_t s;
} partition_param_typdef;

const int16_t icer_wavelet_filter_parameters[][4] = {{0,  4, 4, 0},
                                                     {0,  4, 6, 4},
                                                     {-1, 4, 8, 6},
                                                     {0,  4, 5, 2},
                                                     {0,  3, 8, 6},
                                                     {0,  3, 9, 8},
                                                     {0,  4, 4, 4}};

#define ICER_CONTEXT_MAX 17
#define ICER_DEFAULT_CONTEXT_ZERO_COUNT 2
#define ICER_DEFAULT_CONTEXT_TOTAL_COUNT 4
#define ICER_CONTEXT_RESCALING_CAP 500

enum icer_pixel_contexts {
    ICER_CONTEXT_0 = 0,
    ICER_CONTEXT_1,
    ICER_CONTEXT_2,
    ICER_CONTEXT_3,
    ICER_CONTEXT_4,
    ICER_CONTEXT_5,
    ICER_CONTEXT_6,
    ICER_CONTEXT_7,
    ICER_CONTEXT_8,
    ICER_CONTEXT_9,
    ICER_CONTEXT_10,
    ICER_CONTEXT_11,
    ICER_CONTEXT_12,
    ICER_CONTEXT_13,
    ICER_CONTEXT_14,
    ICER_CONTEXT_15,
    ICER_CONTEXT_16 = ICER_CONTEXT_MAX,
};

const uint8_t icer_context_table_ll_lh_hl[3][3][5] = {
        {
            { ICER_CONTEXT_0,  ICER_CONTEXT_1,  ICER_CONTEXT_2,  ICER_CONTEXT_2,  ICER_CONTEXT_2},
            { ICER_CONTEXT_3,  ICER_CONTEXT_3,  ICER_CONTEXT_3,  ICER_CONTEXT_3,  ICER_CONTEXT_3},
            { ICER_CONTEXT_4,  ICER_CONTEXT_4,  ICER_CONTEXT_4,  ICER_CONTEXT_4,  ICER_CONTEXT_4}
    },
    {
            { ICER_CONTEXT_5,  ICER_CONTEXT_6,  ICER_CONTEXT_7,  ICER_CONTEXT_7,  ICER_CONTEXT_7},
            { ICER_CONTEXT_7,  ICER_CONTEXT_7,  ICER_CONTEXT_7,  ICER_CONTEXT_7,  ICER_CONTEXT_7},
            { ICER_CONTEXT_7,  ICER_CONTEXT_7,  ICER_CONTEXT_7,  ICER_CONTEXT_7,  ICER_CONTEXT_7}
    },
    {
            { ICER_CONTEXT_8,  ICER_CONTEXT_8,  ICER_CONTEXT_8,  ICER_CONTEXT_8,  ICER_CONTEXT_8},
            { ICER_CONTEXT_8,  ICER_CONTEXT_8,  ICER_CONTEXT_8,  ICER_CONTEXT_8,  ICER_CONTEXT_8},
            { ICER_CONTEXT_8,  ICER_CONTEXT_8,  ICER_CONTEXT_8,  ICER_CONTEXT_8,  ICER_CONTEXT_8},
    }
};

const uint8_t icer_context_table_hh[5][5] = {
        { ICER_CONTEXT_0,  ICER_CONTEXT_3,  ICER_CONTEXT_6,  ICER_CONTEXT_8,  ICER_CONTEXT_8},
        { ICER_CONTEXT_1,  ICER_CONTEXT_4,  ICER_CONTEXT_7,  ICER_CONTEXT_8,  ICER_CONTEXT_8},
        { ICER_CONTEXT_2,  ICER_CONTEXT_5,  ICER_CONTEXT_7,  ICER_CONTEXT_8,  ICER_CONTEXT_8},
        { ICER_CONTEXT_2,  ICER_CONTEXT_5,  ICER_CONTEXT_7,  ICER_CONTEXT_8,  ICER_CONTEXT_8},
        { ICER_CONTEXT_2,  ICER_CONTEXT_5,  ICER_CONTEXT_7,  ICER_CONTEXT_8,  ICER_CONTEXT_8},
};

const uint8_t icer_sign_context_table[5][5] = {
        {ICER_CONTEXT_14, ICER_CONTEXT_14, ICER_CONTEXT_15, ICER_CONTEXT_16, ICER_CONTEXT_16},
        {ICER_CONTEXT_14, ICER_CONTEXT_14, ICER_CONTEXT_15, ICER_CONTEXT_16, ICER_CONTEXT_16},
        {ICER_CONTEXT_13, ICER_CONTEXT_13, ICER_CONTEXT_12, ICER_CONTEXT_13, ICER_CONTEXT_13},
        {ICER_CONTEXT_16, ICER_CONTEXT_16, ICER_CONTEXT_15, ICER_CONTEXT_14, ICER_CONTEXT_14},
        {ICER_CONTEXT_16, ICER_CONTEXT_16, ICER_CONTEXT_15, ICER_CONTEXT_14, ICER_CONTEXT_14},
};

/* 1 is negative, 0 is positive */
const uint8_t icer_sign_prediction_table[5][5] = {
        {1, 1, 1, 1, 1},
        {1, 1, 1, 1, 1},
        {0, 0, 0, 1, 1},
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0}
};

enum icer_pixel_categories {
    ICER_CATEGORY_0 = 0,
    ICER_CATEGORY_1,
    ICER_CATEGORY_2,
    ICER_CATEGORY_3 = 3
};

enum icer_subband_types {
    ICER_SUBBAND_LL = 0,
    ICER_SUBBAND_HL,
    ICER_SUBBAND_LH,
    ICER_SUBBAND_HH
};

typedef struct {
    enum icer_subband_types subband_type;
    int16_t zero_count[ICER_CONTEXT_MAX + 1];
    int16_t total_count[ICER_CONTEXT_MAX + 1];
} icer_context_model_typedef;

int icer_wavelet_transform_stages(uint8_t *image, size_t image_w, size_t image_h, uint8_t stages, enum icer_filter_types filt);
int icer_inverse_wavelet_transform_stages(uint8_t *image, size_t image_w, size_t image_h, uint8_t stages, enum icer_filter_types filt);

int icer_wavelet_transform_2d_uint8(uint8_t *image, size_t image_w, size_t image_h, size_t rowstride, enum icer_filter_types filt);
int icer_wavelet_transform_1d_uint8(uint8_t *data, size_t N, size_t stride, enum icer_filter_types filt);

int icer_inverse_wavelet_transform_2d_uint8(uint8_t *image, size_t image_w, size_t image_h, size_t rowstride, enum icer_filter_types filt);
int icer_inverse_wavelet_transform_1d_uint8(uint8_t *data, size_t N, size_t stride, enum icer_filter_types filt);

void icer_reverse_uint8(uint8_t *data, size_t start, size_t end, size_t stride);

void icer_interleave_uint8(uint8_t *data, size_t len, size_t stride);
void icer_deinterleave_uint8(uint8_t *data, size_t len, size_t stride);

uint8_t icer_find_k(size_t len);


void compress_partition_uint8(uint8_t *data, partition_param_typdef *params, size_t rowstride,
                              enum icer_subband_types subband_type, uint8_t lsb);
void compress_bitplane_uint8(uint8_t *data, size_t plane_w, size_t plane_h, size_t rowstride,
                             icer_context_model_typedef *context_model, uint8_t lsb);
void init_context_model_vals(icer_context_model_typedef* context_model, enum icer_subband_types subband_type);
int icer_generate_partition_parameters(partition_param_typdef *params, size_t ll_w, size_t ll_h, uint16_t segments);

static inline uint8_t get_bit_category_uint8(const uint8_t* data, uint8_t lsb);
static inline uint8_t get_bit_category_uint16(const uint16_t* data, uint8_t lsb);
static inline bool get_bit_significance_uint8(const uint8_t* data, uint8_t lsb);
static inline bool get_bit_significance_uint16(const uint16_t* data, uint8_t lsb);
static inline int8_t get_sign_uint8(const uint8_t* data, uint8_t lsb);
static inline int8_t get_sign_uint16(const uint16_t* data, uint8_t lsb);

void icer_to_sign_magnitude_int8(uint8_t *data, size_t len);
void icer_from_sign_magnitude_int8(uint8_t *data, size_t len);

static inline unsigned icer_pow_uint(unsigned base, unsigned exp);

static inline int16_t icer_floor_div_int16(int16_t a, int16_t b);
static inline size_t icer_floor_div_size_t(size_t a, size_t b);
static inline int16_t icer_ceil_div_int16(int16_t a, int16_t b);
static inline size_t icer_ceil_div_size_t(size_t a, size_t b);

static inline int icer_max_int(int a, int b);
static inline unsigned icer_max_uint(unsigned a, unsigned b);

static inline int icer_min_int(int a, int b);


size_t slice_lengths[MAX_K] = {0};

int icer_wavelet_transform_stages(uint8_t *image, size_t image_w, size_t image_h, uint8_t stages,
                                  enum icer_filter_types filt) {
    bool overflow = false;
    size_t low_w = image_w;
    size_t low_h = image_h;

    size_t smallest_w = icer_ceil_div_size_t(image_w, icer_pow_uint(2, stages));
    size_t smallest_h = icer_ceil_div_size_t(image_h, icer_pow_uint(2, stages));

    if (smallest_w < 3 || smallest_h < 3) {
        return ICER_TOO_MANY_STAGES;
    }

    for (uint8_t it = 0; it < stages; it++) {
        overflow |= icer_wavelet_transform_2d_uint8(image, low_w, low_h, image_w, filt);
        low_w = low_w / 2 + low_w % 2;
        low_h = low_h / 2 + low_h % 2;
    }

    icer_to_sign_magnitude_int8(image, image_w * image_h);
    return overflow ? ICER_INTEGER_OVERFLOW : ICER_RESULT_OK;
}

int icer_inverse_wavelet_transform_stages(uint8_t *image, size_t image_w, size_t image_h, uint8_t stages,
                                          enum icer_filter_types filt) {
    bool overflow = false;
    size_t low_w;
    size_t low_h;
    uint8_t decomps;

    size_t smallest_w = icer_ceil_div_size_t(image_w, icer_pow_uint(2, stages));
    size_t smallest_h = icer_ceil_div_size_t(image_h, icer_pow_uint(2, stages));

    if (smallest_w < 3 || smallest_h < 3) {
        return ICER_TOO_MANY_STAGES;
    }

    icer_from_sign_magnitude_int8(image, image_w * image_h);
    for (uint8_t it = 1; it <= stages; it++) {
        decomps = stages - it;
        low_w = icer_ceil_div_size_t(image_w, icer_pow_uint(2, decomps));
        low_h = icer_ceil_div_size_t(image_h, icer_pow_uint(2, decomps));
        overflow |= icer_inverse_wavelet_transform_2d_uint8(image, low_w, low_h, image_w, filt);
    }
    return overflow ? ICER_INTEGER_OVERFLOW : ICER_RESULT_OK;
}

int icer_wavelet_transform_2d_uint8(uint8_t *image, size_t image_w, size_t image_h, size_t rowstride,
                                    enum icer_filter_types filt) {
    bool overflow = false;
    uint8_t *rowstart = image;
    for (size_t r = 0; r < image_h; r++) {
        overflow |= icer_wavelet_transform_1d_uint8(rowstart, image_w, 1, filt);
        rowstart += rowstride;
    }

    uint8_t *colstart = image;
    for (size_t c = 0; c < image_w; c++) {
        overflow |= icer_wavelet_transform_1d_uint8(colstart, image_h, rowstride, filt);
        colstart += 1;
    }

    return overflow ? ICER_INTEGER_OVERFLOW : ICER_RESULT_OK;
}

int icer_inverse_wavelet_transform_2d_uint8(uint8_t *image, size_t image_w, size_t image_h, size_t rowstride,
                                            enum icer_filter_types filt) {
    bool overflow = false;
    uint8_t *colstart = image;
    for (size_t c = 0; c < image_w; c++) {
        overflow |= icer_inverse_wavelet_transform_1d_uint8(colstart, image_h, rowstride, filt);
        colstart += 1;
    }

    uint8_t *rowstart = image;
    for (size_t r = 0; r < image_h; r++) {
        overflow |= icer_inverse_wavelet_transform_1d_uint8(rowstart, image_w, 1, filt);
        rowstart += rowstride;
    }

    return overflow ? ICER_INTEGER_OVERFLOW : ICER_RESULT_OK;
}

static inline int16_t get_r(int8_t *data, size_t n, size_t stride) {
    return (n > 0) ? (int16_t) data[(n - 1) * stride] - (int16_t) data[n * stride] : 1;
}

static inline int16_t get_d(int8_t *data, size_t n, size_t stride, size_t offset, size_t max, bool isodd) {
    return (isodd && n == max) ? 0 : (int16_t) data[(offset + n) * stride];
}

int icer_wavelet_transform_1d_uint8(uint8_t *data, size_t N, size_t stride, enum icer_filter_types filt) {
    size_t low_N, high_N;
    bool is_odd = false;
    bool overflow = false;
    low_N = N / 2 - 1;
    high_N = N / 2 - 1;

    if (N & 1) {
        low_N += 1;
        is_odd = true;
    }

    size_t n1, n2;
    int16_t data1, data2;
    int16_t low, high;

    int8_t *signed_data = (int8_t *) data;
    for (size_t n = 0; n <= low_N; n++) {
        if (!(is_odd && n == low_N)) {
            n1 = (2 * n) * stride;
            n2 = (2 * n + 1) * stride;
            data1 = (int16_t) signed_data[n1];
            data2 = (int16_t) signed_data[n2];

            low = icer_floor_div_int16(data1 + data2, 2);
            high = data1 - data2;

            if (low > INT8_MAX || high > INT8_MAX || low < INT8_MIN || high < INT8_MIN) overflow = true;

            signed_data[n1] = (int8_t) low;
            signed_data[n2] = (int8_t) high;
        } else {
            n1 = (N - 1) * stride;
            data1 = (int16_t) signed_data[n1];
            low = data1;
            high = 0;

            if (low > INT8_MAX || low < INT8_MIN) overflow = true;

            signed_data[n1] = (int8_t) low;
        }
    }

    icer_deinterleave_uint8(data, N, stride);

    int16_t subtract, h;
    size_t offset = low_N + 1;
    bool is_zero = icer_wavelet_filter_parameters[filt][ICER_FILTER_COEF_ALPHA_N1] != 0;
    for (size_t n = 0; n <= high_N; n++) {
        if (n == 0) {
            subtract = icer_floor_div_int16(get_r(signed_data, 1, stride), 4);
        } else if (n == 1 && is_zero) {
            subtract = icer_floor_div_int16(2 * get_r(signed_data, 1, stride)
                                            + 3 * get_r(signed_data, 2, stride)
                                            - 2 * get_d(signed_data, 2, stride, low_N, low_N, is_odd)
                                            + 4, 8);
        } else if (!is_odd && n == N / 2 - 1) {
            subtract = icer_floor_div_int16(get_r(signed_data, N / 2 - 1, stride), 4);
        } else {
            subtract = icer_floor_div_int16(
                    icer_wavelet_filter_parameters[filt][ICER_FILTER_COEF_ALPHA_N1] * get_r(signed_data, n - 1, stride)
                    + icer_wavelet_filter_parameters[filt][ICER_FILTER_COEF_ALPHA_0] * get_r(signed_data, n, stride)
                    + icer_wavelet_filter_parameters[filt][ICER_FILTER_COEF_ALPHA_1] * get_r(signed_data, n + 1, stride)
                    - icer_wavelet_filter_parameters[filt][ICER_FILTER_COEF_BETA] *
                      get_d(signed_data, n + 1, stride, offset, low_N, is_odd)
                    + 8, ICER_FILTER_DENOMINATOR);
        }
        n1 = (offset + n) * stride;
        h = get_d(signed_data, n, stride, offset, low_N, is_odd) - subtract;
        if (h > INT8_MAX || h < INT8_MIN) overflow = true;
        signed_data[n1] = (int8_t) h;
    }

    return overflow ? ICER_INTEGER_OVERFLOW : ICER_RESULT_OK;
}

int icer_inverse_wavelet_transform_1d_uint8(uint8_t *data, size_t N, size_t stride, enum icer_filter_types filt) {
    size_t low_N, high_N;
    bool is_odd = false;
    bool overflow = false;
    low_N = N / 2 - 1;
    high_N = N / 2 - 1;

    if (N & 1) {
        low_N += 1;
        is_odd = true;
    }
    int8_t *signed_data = (int8_t *) data;

    size_t n1;
    int16_t add, d;
    size_t offset = low_N + 1;
    bool is_zero = icer_wavelet_filter_parameters[filt][ICER_FILTER_COEF_ALPHA_N1] != 0;
    for (size_t n, it = 0; it <= high_N; it++) {
        n = high_N - it;
        if (n == 0) {
            add = icer_floor_div_int16(get_r(signed_data, 1, stride), 4);
        } else if (n == 1 && is_zero) {
            add = icer_floor_div_int16(2 * get_r(signed_data, 1, stride)
                                       + 3 * get_r(signed_data, 2, stride)
                                       - 2 * get_d(signed_data, 2, stride, low_N, low_N, is_odd)
                                       + 4, 8);
        } else if (!is_odd && n == N / 2 - 1) {
            add = icer_floor_div_int16(get_r(signed_data, N / 2 - 1, stride), 4);
        } else {
            add = icer_floor_div_int16(
                    icer_wavelet_filter_parameters[filt][ICER_FILTER_COEF_ALPHA_N1] * get_r(signed_data, n - 1, stride)
                    + icer_wavelet_filter_parameters[filt][ICER_FILTER_COEF_ALPHA_0] * get_r(signed_data, n, stride)
                    + icer_wavelet_filter_parameters[filt][ICER_FILTER_COEF_ALPHA_1] * get_r(signed_data, n + 1, stride)
                    - icer_wavelet_filter_parameters[filt][ICER_FILTER_COEF_BETA] *
                      get_d(signed_data, n + 1, stride, offset, low_N, is_odd)
                    + 8, ICER_FILTER_DENOMINATOR);
        }
        n1 = (offset + n) * stride;
        d = get_d(signed_data, n, stride, offset, low_N, is_odd) + add;
        if (d > INT8_MAX || d < INT8_MIN) overflow = true;

        signed_data[n1] = (int8_t) d;
    }

    size_t n2;
    int16_t low, high, tmp;
    for (size_t n = 0; n <= low_N; n++) {
        if (!(is_odd && n == low_N)) {
            n1 = n * stride;
            n2 = (offset + n) * stride;

            low = (int16_t) signed_data[n1];
            high = (int16_t) signed_data[n2];

            tmp = low + icer_floor_div_int16(high + 1, 2);

            if (tmp > INT8_MAX || (tmp - high) > INT8_MAX || tmp < INT8_MIN || (tmp - high) < INT8_MIN) overflow = true;

            signed_data[n1] = (int8_t) tmp;
            signed_data[n2] = (int8_t) (tmp - high);
        } else {
            n1 = n * stride;

            low = (int16_t) signed_data[n1];
            high = 0;

            tmp = low + icer_floor_div_int16(high + 1, 2);

            if (tmp > INT8_MAX || tmp < INT8_MIN) overflow = true;

            signed_data[n1] = (int8_t) tmp;
        }
    }

    icer_interleave_uint8(data, N, stride);

    return overflow ? ICER_INTEGER_OVERFLOW : ICER_RESULT_OK;
}

void icer_reverse_uint8(uint8_t *data, size_t start, size_t end, size_t stride) {
    size_t n1, n2;
    uint8_t swap;
    while (start < end) {
        n1 = start * stride;
        n2 = end * stride;
        swap = data[n1];
        data[n1] = data[n2];
        data[n2] = swap;

        start++;
        end--;
    }
}

void icer_interleave_uint8(uint8_t *data, size_t len, size_t stride) {
    size_t n = len;
    size_t processed = 0, left, halfleft;
    size_t segment, halfseg;
    uint8_t k;

    size_t j, n1;
    uint8_t swap, num;

    bool is_odd = len & 1;

    while (processed < n) {
        k = icer_find_k(n - processed);
        segment = slice_lengths[k];
        halfseg = segment / 2;

        left = n - processed;
        halfleft = left / 2 - (is_odd ? 0 : 1);

        if (left > 0) {
            icer_reverse_uint8(data, processed + halfseg, processed + halfleft + halfseg, stride);
            icer_reverse_uint8(data, processed + halfseg, processed + segment - 1, stride);
            icer_reverse_uint8(data, processed + segment, processed + halfleft + halfseg, stride);
        }

        for (size_t i = 1; i < segment; i *= 3) {
            j = i;
            n1 = (j + processed) * stride;
            num = data[n1];
            do {
                if (j < halfseg) {
                    j *= 2;
                } else {
                    j = (j - halfseg) * 2 + 1;
                }

                n1 = (j + processed) * stride;

                swap = data[n1];
                data[n1] = num;
                num = swap;
            } while (j != i);
        }

        processed += segment;
    }
}

void icer_deinterleave_uint8(uint8_t *data, size_t len, size_t stride) {
    size_t n = len;
    size_t processed = 0;
    size_t segment, halfseg;
    uint8_t k;

    size_t j, n1;
    uint8_t swap, num;

    while (processed < n) {
        k = icer_find_k(n - processed);
        segment = slice_lengths[k];
        halfseg = segment / 2;

        for (size_t i = 1; i < segment; i *= 3) {
            j = i;
            n1 = (j + processed) * stride;
            num = data[n1];
            do {
                if (j & 1) {
                    j = halfseg + j / 2;
                } else {
                    j /= 2;
                }

                n1 = (j + processed) * stride;

                swap = data[n1];
                data[n1] = num;
                num = swap;
            } while (j != i);
        }

        if (processed > 0) {
            icer_reverse_uint8(data, processed / 2, processed - 1, stride);
            icer_reverse_uint8(data, processed, processed + halfseg - 1, stride);
            icer_reverse_uint8(data, processed / 2, processed + halfseg - 1, stride);
        }

        processed += segment;
    }
}

uint8_t icer_find_k(size_t len) {
    uint8_t max_k = MAX_K - 1;
    uint8_t min_k = 0;

    if (slice_lengths[min_k] == 0) slice_lengths[min_k] = icer_pow_uint(3, min_k) + 1;
    if (slice_lengths[max_k] == 0) slice_lengths[max_k] = icer_pow_uint(3, max_k) + 1;

    uint8_t mid;
    uint8_t res = 0;
    while (min_k < max_k) {
        mid = (max_k + min_k) / 2;
        if (slice_lengths[mid] == 0) slice_lengths[mid] = icer_pow_uint(3, mid) + 1;

        if (len > slice_lengths[mid]) {
            min_k = mid + 1;
            res = mid;
        } else if (len < slice_lengths[mid]) {
            max_k = mid - 1;
        } else {
            break;
        }
    }

    return res;
}

void compress_partition_uint8(uint8_t *data, partition_param_typdef *params, size_t rowstride,
                              enum icer_subband_types subband_type, uint8_t lsb) {
    size_t segment_w, segment_h;
    uint8_t *segment_start;
    uint16_t segment_num = 0;

    size_t partition_col_ind;
    size_t partition_row_ind = 0;

    icer_context_model_typedef context_model;

    /*
     * process top region which consists of c columns
     * height of top region is h_t and it contains r_t rows
     */
    for (uint16_t row = 0; row < params->r_t; row++) {
        /*
         * the first r_t0 rows have height y_t
         * the remainder have height y_t + 1
         */
        segment_h = params->y_t + ((row > params->r_t0) ? 1 : 0);
        partition_col_ind = 0;

        for (uint16_t col = 0; col < params->c; col++) {
            /* the first c_t0 columns have width x_t
             * the remainder have width x_t + 1
             */
            segment_w = params->x_t + ((col > params->c_t0) ? 1 : 0);
            segment_start = data + partition_row_ind * rowstride + partition_col_ind;
            segment_num++;
            partition_col_ind += segment_w;

            init_context_model_vals(&context_model, subband_type);
            printf("t segment no: %d\n", segment_num);
            compress_bitplane_uint8(segment_start, segment_w, segment_h, rowstride, &context_model, lsb);
        }
        partition_row_ind += segment_h;
    }

    /*
     * if the bottom region exists, process bottom region
     * which consists of c+1 columns
     */
    for (uint16_t row = 0; row < (params->r - params->r_t); row++) {
        /*
         * the first r_b0 rows have height y_b
         * the remainder have height y_b + 1
         */
        segment_h = params->y_b + ((row > params->r_b0) ? 1 : 0);
        partition_col_ind = 0;

        for (uint16_t col = 0; col < (params->c + 1); col++) {
            /* the first c_b0 columns have width x_b
             * the remainder have width x_b + 1
             */
            segment_w = params->x_b + ((col > params->c_b0) ? 1 : 0);
            segment_start = data + partition_row_ind * rowstride + partition_col_ind;
            segment_num++;
            partition_col_ind += segment_w;

            init_context_model_vals(&context_model, subband_type);
            printf("b segment no: %d\n", segment_num);
            compress_bitplane_uint8(segment_start, segment_w, segment_h, rowstride, &context_model, lsb);
        }
        partition_row_ind += segment_h;
    }
}

void enc_bit(uint8_t bit, int16_t zeros, int16_t total) {
    printf("%1u %f\n", bit & 1, (float)zeros/(float)total);
}

void compress_bitplane_uint8(uint8_t *data, size_t plane_w, size_t plane_h, size_t rowstride,
                             icer_context_model_typedef *context_model, uint8_t lsb) {
    uint8_t *pos;
    uint8_t *rowstart = data;
    int category;
    bool bit;
    uint8_t mask = 0b1 << (lsb-1);

    uint8_t *h0, *h1, *v0, *v1, *d0, *d1, *d2, *d3;
    uint8_t h, v, d, tmp;
    int8_t sh0, sh1, sv0, sv1;
    uint8_t sh, sv;
    uint8_t pred_sign, actual_sign, agreement_bit;
    enum icer_pixel_contexts sign_context;
    size_t vert_bound = plane_h - 1;
    size_t hor_bound = plane_w - 1;

    uint8_t prev_plane = lsb+1;

    for (size_t row = 0; row < plane_h; row++) {
        pos = rowstart;
        /* keep track of 8 adjacent pixels */
        h0 = pos - 1; h1 = pos + 1;
        v0 = pos - rowstride; v1 = pos + rowstride;
        d0 = v0 - 1; d1 = v0 + 1;
        d2 = v1 - 1; d3 = v1 + 1;

        enum icer_pixel_contexts cntxt;

        for (size_t col = 0; col < plane_w; col++) {
            category = get_bit_category_uint8(pos, lsb);
            bit = ((*pos) & mask) != 0;

            if (category == ICER_CATEGORY_3) {
                /* pass to uncoded bin */

            } else {
                if (category == ICER_CATEGORY_0 || category == ICER_CATEGORY_1) {
                    h = 0;
                    v = 0;
                    d = 0;

                    /* consider the horizontally adjacent pixels */
                    if (col > 0) h += get_bit_significance_uint8(h0, lsb) & 1;
                    if (col < hor_bound) h += get_bit_significance_uint8(h1, prev_plane) & 1;

                    /* consider the vertically adjacent pixels */
                    if (row > 0) v += get_bit_significance_uint8(v0, lsb) & 1;
                    if (row < vert_bound) v += get_bit_significance_uint8(v1, prev_plane) & 1;

                    /* consider the diagonally adjacent pixels */
                    if (col > 0 && row > 0) d += get_bit_significance_uint8(d0, lsb) & 1;
                    if (col > 0 && row < vert_bound) d += get_bit_significance_uint8(d2, prev_plane) & 1;
                    if (col < hor_bound && row > 0) d += get_bit_significance_uint8(d1, lsb) & 1;
                    if (col < hor_bound && row < vert_bound) d += get_bit_significance_uint8(d3, prev_plane) & 1;
                }

                if (category == ICER_CATEGORY_0) {
                    if (context_model->subband_type == ICER_SUBBAND_HL) {
                        tmp = h;
                        h = v;
                        v = tmp;
                    }

                    if (context_model->subband_type != ICER_SUBBAND_HH) {
                        cntxt = icer_context_table_ll_lh_hl[h][v][d];
                    } else {
                        cntxt = icer_context_table_hh[h + v][d];
                    }
                } else if (category == ICER_CATEGORY_1) {
                    cntxt = (h + v == 0) ? ICER_CONTEXT_9 : ICER_CONTEXT_10;
                } else if (category == ICER_CATEGORY_2) {
                    cntxt = ICER_CONTEXT_11;
                }

                //todo: encode bit
                enc_bit(bit, context_model->zero_count[cntxt], context_model->total_count[cntxt]);

                context_model->total_count[cntxt]++;
                context_model->zero_count[cntxt] += (int16_t)(!bit);
                if (context_model->total_count[cntxt] >= ICER_CONTEXT_RESCALING_CAP) {
                    context_model->total_count[cntxt] /= 2;
                    context_model->zero_count[cntxt] /= 2;

                    //todo: when necessary, the count of zeros is rounded in the direction that makes the probability-of-zero estimate closer to 1/2)
                }

                if (category == ICER_CATEGORY_0 && bit) {
                    /* bit is the first magnitude bit to be encoded, thus we will encode the sign bit */
                    /* predict sign bit */
                    sh0 = 0; sh1 = 0;
                    sv0 = 0; sv1 = 0;

                    /* consider the horizontally adjacent pixels */
                    //if (col > 0) sh0 = get_sign_uint8(h0, lsb);
                    //if (col < hor_bound) sh1 = get_sign_uint8(h1, prev_plane);

                    /* consider the vertically adjacent pixels */
                    //if (row > 0) sv0 = get_sign_uint8(v0, lsb) & 1;
                    //if (row < vert_bound) sv1 = get_sign_uint8(v1, prev_plane) & 1;

                    sh = sh0 + sh1 + 2;
                    sv = sv0 + sv1 + 2;

                    if (context_model->subband_type == ICER_SUBBAND_HL) {
                        tmp = sh;
                        sh = sv;
                        sv = tmp;
                    }

                    sign_context = icer_sign_context_table[sh][sv];
                    pred_sign = icer_sign_prediction_table[sh][sv];
                    actual_sign = (*pos) & 0x80;

                    agreement_bit = (pred_sign ^ actual_sign) & 1;

                    //todo: encode sign bit
                    enc_bit(agreement_bit, context_model->zero_count[sign_context], context_model->total_count[sign_context]);

                    context_model->total_count[sign_context]++;
                    context_model->zero_count[sign_context] += (int16_t)(agreement_bit == 0);
                    if (context_model->total_count[sign_context] >= ICER_CONTEXT_RESCALING_CAP) {
                        context_model->total_count[sign_context] /= 2;
                        context_model->zero_count[sign_context] /= 2;

                        //todo: when necessary, the count of zeros is rounded in the direction that makes the probability-of-zero estimate closer to 1/2)
                    }
                }
            }

            pos++;
            /* increment position pointers of 8 adjacent pixels */
            h0++; h1++; v0++; v1++; d0++; d1++; d2++; d2++;
        }
        rowstart += rowstride;
    }
}

void init_context_model_vals(icer_context_model_typedef* context_model, enum icer_subband_types subband_type) {
    context_model->subband_type = subband_type;
    for (size_t i = 0;i <= ICER_CONTEXT_MAX;i++) {
        context_model->zero_count[i] = ICER_DEFAULT_CONTEXT_ZERO_COUNT;
        context_model->total_count[i] = ICER_DEFAULT_CONTEXT_TOTAL_COUNT;
    }
}

int icer_generate_partition_parameters(partition_param_typdef *params, size_t ll_w, size_t ll_h, uint16_t segments) {
    uint16_t r, c, r_t, h_t, x_t, c_t0, y_t, r_t0;
    uint16_t x_b = 0, c_b0 = 0, y_b = 0, r_b0 = 0;

    if (segments > (ll_w * ll_h) || segments > ICER_MAX_SEGMENTS) {
        return ICER_TOO_MANY_SEGMENTS;
    }

    if (ll_h > (segments - 1) * ll_w) {
        r = segments;
    } else {
        for (r = 1; r < segments && (r + 1) * r * ll_w < ll_h * segments; r++);
    }
    c = segments / r;
    r_t = (c + 1) * r - segments;
    h_t = icer_max_uint(r_t, ((2 * ll_h * c * r_t + segments) / 2) / segments);
    x_t = ll_w / c;
    c_t0 = (x_t + 1) * c - ll_w;
    y_t = h_t / r_t;
    r_t0 = (y_t + 1) * r_t - h_t;

    if (r_t < r) {
        x_b = ll_w / (c + 1);
        c_b0 = (x_b + 1) * (c + 1) - ll_w;
        y_b = (ll_h - h_t) / (r - r_t);
        r_b0 = (y_b + 1) * (r - r_t) - (ll_h - h_t);
    }

    params->w = ll_w;
    params->h = ll_h;
    params->s = segments;

    params->r = r;
    params->c = c;
    params->r_t = r_t;
    params->h_t = h_t;
    params->x_t = x_t;
    params->c_t0 = c_t0;
    params->y_t = y_t;
    params->r_t0 = r_t0;

    params->x_b = x_b;
    params->c_b0 = c_b0;
    params->y_b = y_b;
    params->r_b0 = r_b0;

    return ICER_RESULT_OK;
}

static inline uint8_t get_bit_category_uint8(const uint8_t* data, uint8_t lsb) {
    return icer_min_int(__builtin_popcount(((*data) & 0x7f) >> lsb), 3);
}

static inline uint8_t get_bit_category_uint16(const uint16_t* data, uint8_t lsb) {
    return icer_min_int(__builtin_popcount(((*data) & 0x7fff) >> lsb), 3);
}

static inline bool get_bit_significance_uint8(const uint8_t* data, uint8_t lsb) {
    return __builtin_popcount(((*data) & 0x7f) >> lsb) != 0;
}

static inline bool get_bit_significance_uint16(const uint16_t* data, uint8_t lsb) {
    return __builtin_popcount(((*data) & 0x7fff) >> lsb) != 0;
}

static inline int8_t get_sign_uint8(const uint8_t* data, uint8_t lsb) {
    return ((int8_t)(*data) >> 7) * (int8_t)get_bit_significance_uint8(data, lsb);
}

static inline int8_t get_sign_uint16(const uint16_t* data, uint8_t lsb) {
    return ((int16_t)(*data) >> 15) * (int8_t)get_bit_significance_uint16(data, lsb);
}

/* this assumes that the platform you are running the code on stores numbers as 2's complement */
void icer_to_sign_magnitude_int8(uint8_t *data, size_t len) {
    uint8_t mask;
    for (uint8_t *it = data; it < data + len; it++) {
        mask = (int8_t) (*it) >> 7;
        (*it) = (((int8_t) (*it) + (int8_t) mask) ^ mask) | ((*it) & 0b10000000);
    }
}

/* this assumes that the platform you are running the code on stores numbers as 2's complement */
void icer_from_sign_magnitude_int8(uint8_t *data, size_t len) {
    uint8_t mask;
    for (uint8_t *it = data; it < data + len; it++) {
        mask = (int8_t) (*it) >> 7;
        (*it) = (~mask & (*it)) | (((int8_t) ((*it) & 0b10000000) - (int8_t) (*it)) & mask);
    }
}

static inline unsigned icer_pow_uint(unsigned base, unsigned exp) {
    unsigned res = 1;
    while (exp) {
        if (exp & 1) res *= base;
        exp /= 2;
        base *= base;
    }
    return res;
}

static inline int16_t icer_floor_div_int16(int16_t a, int16_t b) {
    int16_t d = a / b;
    int16_t r = a %
                b; /* optimizes into single division, as per https://stackoverflow.com/questions/46265403/fast-floor-of-a-signed-integer-division-in-c-c */
    return r ? (d - ((a < 0) ^ (b < 0))) : d;
}

static inline size_t icer_floor_div_size_t(size_t a, size_t b) {
    size_t d = a / b;
    size_t r = a %
                b; /* optimizes into single division, as per https://stackoverflow.com/questions/46265403/fast-floor-of-a-signed-integer-division-in-c-c */
    return r ? (d - ((a < 0) ^ (b < 0))) : d;
}

static inline int16_t icer_ceil_div_int16(int16_t a, int16_t b) {
    int16_t d = a / b;
    int16_t r = a %
                b; /* optimizes into single division, as per https://stackoverflow.com/questions/46265403/fast-floor-of-a-signed-integer-division-in-c-c */
    return r ? d + 1 : d;
}

static inline size_t icer_ceil_div_size_t(size_t a, size_t b) {
    size_t d = a / b;
    size_t r = a %
               b; /* optimizes into single division, as per https://stackoverflow.com/questions/46265403/fast-floor-of-a-signed-integer-division-in-c-c */
    return r ? d + 1 : d;
}

static inline int icer_max_int(int a, int b) {
    return (a > b) ? a : b;
}

static inline unsigned icer_max_uint(unsigned a, unsigned b) {
    return (a > b) ? a : b;
}

static inline int icer_min_int(int a, int b) {
    return (a < b) ? a : b;
}

#endif //ICER_COMPRESSION_ICER_H
