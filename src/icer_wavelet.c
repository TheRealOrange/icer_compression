//
// Created by linyi on 14/3/2023.
//

#include "icer.h"

int icer_wavelet_transform_stages(uint8_t *image, size_t image_w, size_t image_h, uint8_t stages,
                                  enum icer_filter_types filt) {
    bool overflow = false;
    size_t low_w = image_w;
    size_t low_h = image_h;

    size_t smallest_w = get_dim_n_low_stages(image_w, stages);
    size_t smallest_h = get_dim_n_low_stages(image_h, stages);

    if (smallest_w < 3 || smallest_h < 3) {
        return ICER_TOO_MANY_STAGES;
    }

    for (uint8_t it = 0; it < stages; it++) {
        overflow |= icer_wavelet_transform_2d_uint8(image, low_w, low_h, image_w, filt);
        low_w = low_w / 2 + low_w % 2;
        low_h = low_h / 2 + low_h % 2;
    }

    return overflow ? ICER_INTEGER_OVERFLOW : ICER_RESULT_OK;
}

int icer_inverse_wavelet_transform_stages(uint8_t *image, size_t image_w, size_t image_h, uint8_t stages,
                                          enum icer_filter_types filt) {
    bool overflow = false;
    size_t low_w;
    size_t low_h;
    uint8_t decomps;

    size_t smallest_w = get_dim_n_low_stages(image_w, stages);
    size_t smallest_h = get_dim_n_low_stages(image_h, stages);

    if (smallest_w < 3 || smallest_h < 3) {
        return ICER_TOO_MANY_STAGES;
    }

    //icer_from_sign_magnitude_int8(image, image_w * image_h);
    for (uint8_t it = 1; it <= stages; it++) {
        decomps = stages - it;
        low_w = get_dim_n_low_stages(image_w, decomps);
        low_h = get_dim_n_low_stages(image_h, decomps);
        overflow |= icer_inverse_wavelet_transform_2d_uint8(image, low_w, low_h, image_w, filt);
    }
    return overflow ? ICER_INTEGER_OVERFLOW : ICER_RESULT_OK;
}

size_t get_dim_n_low_stages(size_t dim, uint8_t stages) {
    return icer_ceil_div_size_t(dim, icer_pow_uint(2, stages));
}

size_t get_dim_n_high_stages(size_t dim, uint8_t stages) {
    return icer_floor_div_size_t(dim, icer_pow_uint(2, stages));
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

/* interleaving and deinterleaving helper functions */

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