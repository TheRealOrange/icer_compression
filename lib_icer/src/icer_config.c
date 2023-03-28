//
// Created by linyi on 19/3/2023.
//
#include "icer.h"

#ifdef USE_ENCODE_FUNCTIONS
icer_custom_code_typedef icer_custom_coding_scheme[ICER_ENCODER_BIN_MAX + 1][CUSTOM_CODING_MAX_LOOKUP];
#endif

#ifdef USE_DECODE_FUNCTIONS
icer_custom_code_typedef icer_custom_decode_scheme[ICER_ENCODER_BIN_MAX + 1][CUSTOM_CODING_MAX_LOOKUP];
#endif

icer_custom_flush_typedef icer_custom_code_flush_bits[ICER_ENCODER_BIN_MAX + 1][CUSTOM_CODE_FLUSH_MAX_LOOKUP + 1][MAX_NUM_BITS_BEFORE_FLUSH + 1];

icer_golomb_code_typedef icer_golomb_coders[ICER_ENCODER_BIN_MAX + 1];

const int16_t icer_wavelet_filter_parameters[][4] = {{0,  4, 4, 0},
                                                     {0,  4, 6, 4},
                                                     {-1, 4, 8, 6},
                                                     {0,  4, 5, 2},
                                                     {0,  3, 8, 6},
                                                     {0,  3, 9, 8},
                                                     {0,  4, 4, 4}};

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

const uint32_t icer_bin_probability_cutoffs[ICER_ENCODER_BIN_MAX+1] = {
        35298,
        37345,
        40503,
        43591,
        47480,
        50133,
        53645,
        55902,
        57755,
        58894,
        60437,
        62267,
        63613,
        64557,
        65134,
        65392,
        65536
};

const int32_t icer_bin_coding_scheme[ICER_ENCODER_BIN_MAX+1] = {
        0,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        5,
        6,
        7,
        11,
        17,
        31,
        70,
        200,
        512
};
size_t icer_slice_lengths[MAX_K] = {0};
