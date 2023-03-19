#ifndef ICER_COMPRESSION_ICER_H
#define ICER_COMPRESSION_ICER_H

#include <stdbool.h>
#include <stdlib.h>

#include "crc.h"

#ifdef _MSC_VER
#include <intrin.h>

#define __builtin_popcount __popcnt
#define __builtin_clz __clz

uint32_t __inline __clz( uint32_t value ) {
    DWORD leading_zero = 0;

    if ( _BitScanReverse( &leading_zero, value ) ) {
       return 31 - leading_zero;
    } else {
         // Same remarks as above
         return 32;
    }
}
#endif


#define MAX_K 12

enum icer_status {
    ICER_RESULT_OK = 0,
    ICER_INTEGER_OVERFLOW = 1,
    ICER_SIZE_ERROR = 2,
    ICER_TOO_MANY_SEGMENTS = 3,
    ICER_TOO_MANY_STAGES = 4,
    ICER_BYTE_QUOTA_EXCEEDED = 5,
    ICER_BITPLANE_OUT_OF_RANGE
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

#define ICER_CONTEXT_MAX 16
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
    uint32_t zero_count[ICER_CONTEXT_MAX + 1];
    uint32_t total_count[ICER_CONTEXT_MAX + 1];
} icer_context_model_typedef;

#define ICER_ENCODER_BIN_MAX 16

enum icer_encoder_bins {
    ICER_ENC_BIN_1 = 0,
    ICER_ENC_BIN_2,
    ICER_ENC_BIN_3,
    ICER_ENC_BIN_4,
    ICER_ENC_BIN_5,
    ICER_ENC_BIN_6,
    ICER_ENC_BIN_7,
    ICER_ENC_BIN_8,
    ICER_ENC_BIN_9,
    ICER_ENC_BIN_10,
    ICER_ENC_BIN_11,
    ICER_ENC_BIN_12,
    ICER_ENC_BIN_13,
    ICER_ENC_BIN_14,
    ICER_ENC_BIN_15,
    ICER_ENC_BIN_16,
    ICER_ENC_BIN_17 = ICER_ENCODER_BIN_MAX
};

#define ICER_ENC_BUF_BITS_OFFSET 11
#define ICER_ENC_BUF_SHIFTOUT_OFFSET 6
#define ICER_ENC_BUF_DONE_MASK 0b0000010000000000
#define ICER_ENC_BUF_DATA_MASK 0b0000001111111111
#define ICER_ENC_BUF_INFO_MASK 0b1111100000000000

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

#define ICER_BIN_PROBABILITY_DENOMINATOR 65536

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

typedef struct {
    uint8_t input_code_bits;
    uint8_t output_code_bits;
    uint8_t output_code;
} custom_code_typedef;

typedef struct {
    uint8_t flush_bit;
    uint8_t flush_bit_numbers;
} custom_flush_typedef;

#define CUSTOM_CODING_MAX_LOOKUP 32
custom_code_typedef custom_coding_scheme[ICER_ENCODER_BIN_MAX+1][CUSTOM_CODING_MAX_LOOKUP];

#define CUSTOM_CODE_FLUSH_MAX_LOOKUP 8
#define MAX_NUM_BITS_BEFORE_FLUSH    5
custom_flush_typedef custom_code_flush_bits[ICER_ENCODER_BIN_MAX+1][CUSTOM_CODE_FLUSH_MAX_LOOKUP+1][MAX_NUM_BITS_BEFORE_FLUSH+1];

typedef struct {
    uint16_t m;
    uint16_t l;
    uint16_t i;
} golomb_code_typedef;

golomb_code_typedef golomb_coders[ICER_ENCODER_BIN_MAX+1];

typedef struct {
    uint16_t subband_number;
    uint8_t subband_type;
    uint8_t ll_mean_val;
    uint8_t lsb;
} packet_context;

#define ICER_PACKET_PREAMBLE 0x605B

typedef struct {
    uint16_t preamble;
    uint16_t subband_number;
    uint8_t subband_type;
    uint8_t segment_number;
    uint8_t lsb;
    uint8_t ll_mean_val;
    uint32_t data_length; // store data length in bits for the decoder
    uint32_t crc32;
    uint8_t *data;
} image_segment_typedef;

typedef struct {
    size_t size_used;
    size_t size_allocated;
    uint8_t *data_start;
} output_data_buf_typedef;

typedef struct {
    size_t buffer_length;
    size_t head;
    size_t tail;
    size_t used;
    size_t output_ind;
    size_t max_output_length;
    uint16_t *encode_buffer;
    uint8_t output_bit_offset;
    uint8_t *output_buffer;
    int16_t bin_current_buf[ICER_ENCODER_BIN_MAX+1];
    int16_t bin_current_buf_bits[ICER_ENCODER_BIN_MAX+1];
} encoder_context_typedef;

typedef struct {
    size_t decoded_words;
    size_t encode_ind;
    uint8_t bit_offset;
    uint8_t *encoded_words;
    int16_t bin_current_buf[ICER_ENCODER_BIN_MAX+1];
    int16_t bin_current_buf_bits[ICER_ENCODER_BIN_MAX+1];
    int16_t bin_decode_index[ICER_ENCODER_BIN_MAX+1];
} decoder_context_typedef;

int icer_init();

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


int compress_partition_uint8(uint8_t *data, partition_param_typdef *params, size_t rowstride, packet_context *pkt_context,
                         output_data_buf_typedef *output_data);
int compress_bitplane_uint8(uint8_t *data, size_t plane_w, size_t plane_h, size_t rowstride,
                            icer_context_model_typedef *context_model,
                            encoder_context_typedef *encoder_context,
                            packet_context *pkt_context);

int icer_encode_bit(encoder_context_typedef *encoder_context, uint8_t bit, uint32_t zero_cnt, uint32_t total_cnt);
int icer_decode_bit(decoder_context_typedef *decoder_context, uint8_t *bit, uint32_t zero_cnt, uint32_t total_cnt);

int flush_encode(encoder_context_typedef *encoder_context);
void init_context_model_vals(icer_context_model_typedef* context_model, enum icer_subband_types subband_type);
void init_entropy_coder_context(encoder_context_typedef *encoder_context, uint16_t *encode_buffer, size_t buffer_length, uint8_t *encoder_out, size_t enc_out_max);
int icer_generate_partition_parameters(partition_param_typdef *params, size_t ll_w, size_t ll_h, uint16_t segments);

uint32_t icer_calculate_packet_crc32(image_segment_typedef *pkt);
int icer_allocate_data_packet(image_segment_typedef **pkt, output_data_buf_typedef *output_data, uint8_t segment_num, packet_context *context);

static inline uint8_t get_bit_category_uint8(const uint8_t* data, uint8_t lsb);
static inline uint8_t get_bit_category_uint16(const uint16_t* data, uint8_t lsb);
static inline bool get_bit_significance_uint8(const uint8_t* data, uint8_t lsb);
static inline bool get_bit_significance_uint16(const uint16_t* data, uint8_t lsb);
static inline int8_t get_sign_uint8(const uint8_t* data, uint8_t lsb);
static inline int8_t get_sign_uint16(const uint16_t* data, uint8_t lsb);

void icer_init_output_struct(output_data_buf_typedef *out, uint8_t *data, size_t len);
static inline uint16_t pop_buf(encoder_context_typedef *cntxt);
static inline int16_t alloc_buf(encoder_context_typedef *cntxt);
int icer_compute_bin(uint32_t zero_cnt, uint32_t total_cnt);

void icer_to_sign_magnitude_int8(uint8_t *data, size_t len);
void icer_from_sign_magnitude_int8(uint8_t *data, size_t len);

static inline unsigned icer_pow_uint(unsigned base, unsigned exp);

static inline int16_t icer_floor_div_int16(int16_t a, int16_t b);
static inline size_t icer_floor_div_size_t(size_t a, size_t b);
static inline int16_t icer_ceil_div_int16(int16_t a, int16_t b);
static inline uint32_t icer_ceil_div_uint32(uint32_t a, uint32_t b);
static inline size_t icer_ceil_div_size_t(size_t a, size_t b);

static inline int icer_max_int(int a, int b);
static inline unsigned icer_max_uint(unsigned a, unsigned b);

static inline int icer_min_int(int a, int b);

size_t slice_lengths[MAX_K] = {0};

#define INIT_CODING_SCHEME(bin, inp, inp_bits, out, out_bits) ({ \
custom_coding_scheme[bin][inp].input_code_bits = inp_bits;       \
custom_coding_scheme[bin][inp].output_code = out;                \
custom_coding_scheme[bin][inp].output_code_bits = out_bits;      \
})

#define INIT_FLUSH_BITS(bin, inp, inp_bits, out, out_bits) ({ \
custom_code_flush_bits[bin][inp][inp_bits].flush_bit = out;   \
custom_code_flush_bits[bin][inp][inp_bits].flush_bit_numbers = out_bits; \
})

int icer_init() {
    for (int it = 0;it <= ICER_ENCODER_BIN_MAX;it++) {
        if (icer_bin_coding_scheme[it] > 0) {
            unsigned int m = icer_bin_coding_scheme[it];
            golomb_coders[it].m = m;

            // compute ceil( log2( m ) )
            unsigned int l = 31-__builtin_clz(m);
            l += ((m ^ (1 << l)) != 0);

            // compute 2^l - m
            unsigned int i = icer_pow_uint(2, l) - m;

            golomb_coders[it].i = i;
            golomb_coders[it].l = l;
        }
    }

    for (int it = 0;it <= ICER_ENCODER_BIN_MAX;it++) {
        for (int j = 0;j < CUSTOM_CODING_MAX_LOOKUP;j++) custom_coding_scheme[it][j].input_code_bits = 0;
    }

    INIT_CODING_SCHEME(ICER_ENC_BIN_2,    0b01, 2,    0b01, 2);
    INIT_CODING_SCHEME(ICER_ENC_BIN_2,   0b011, 3,   0b110, 3);
    INIT_CODING_SCHEME(ICER_ENC_BIN_2,  0b0111, 4,  0b1111, 4);
    INIT_CODING_SCHEME(ICER_ENC_BIN_2,  0b1111, 4, 0b00001, 5);
    INIT_CODING_SCHEME(ICER_ENC_BIN_2,    0b10, 2,    0b10, 2);
    INIT_CODING_SCHEME(ICER_ENC_BIN_2,   0b100, 3,   0b001, 3);
    INIT_CODING_SCHEME(ICER_ENC_BIN_2,  0b1000, 4,  0b0001, 4);
    INIT_CODING_SCHEME(ICER_ENC_BIN_2, 0b00001, 5, 0b00000, 5);
    INIT_CODING_SCHEME(ICER_ENC_BIN_2, 0b00000, 5,  0b1110, 4);

    INIT_FLUSH_BITS(ICER_ENC_BIN_2,      0b1, 1, 0, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_2,     0b11, 2, 0, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_2,    0b111, 3, 0, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_2,      0b0, 1, 1, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_2,     0b00, 2, 1, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_2,    0b000, 3, 1, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_2,   0b0000, 4, 0, 1);

    INIT_CODING_SCHEME(ICER_ENC_BIN_3,    0b10, 2,    0b10, 2);
    INIT_CODING_SCHEME(ICER_ENC_BIN_3,   0b100, 3,    0b00, 2);
    INIT_CODING_SCHEME(ICER_ENC_BIN_3,  0b0000, 4,   0b110, 3);
    INIT_CODING_SCHEME(ICER_ENC_BIN_3, 0b00011, 5, 0b01001, 5);
    INIT_CODING_SCHEME(ICER_ENC_BIN_3, 0b00010, 5,  0b1111, 4);
    INIT_CODING_SCHEME(ICER_ENC_BIN_3,    0b01, 2,   0b011, 3);
    INIT_CODING_SCHEME(ICER_ENC_BIN_3,  0b0011, 4,  0b1110, 4);
    INIT_CODING_SCHEME(ICER_ENC_BIN_3,  0b1011, 4, 0b01000, 5);
    INIT_CODING_SCHEME(ICER_ENC_BIN_3,   0b111, 3,  0b0101, 4);

    INIT_FLUSH_BITS(ICER_ENC_BIN_3,      0b0, 1, 1, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_3,     0b00, 2, 1, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_3,    0b000, 3, 0, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_3,   0b1000, 4, 0, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_3,      0b1, 1, 0, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_3,     0b11, 2, 1, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_3,    0b011, 3, 0, 1);

    INIT_CODING_SCHEME(ICER_ENC_BIN_4,    0b10, 2,    0b01, 2);
    INIT_CODING_SCHEME(ICER_ENC_BIN_4,   0b100, 3,   0b110, 3);
    INIT_CODING_SCHEME(ICER_ENC_BIN_4,   0b000, 3,    0b00, 2);
    INIT_CODING_SCHEME(ICER_ENC_BIN_4,    0b01, 2,    0b10, 2);
    INIT_CODING_SCHEME(ICER_ENC_BIN_4,    0b11, 2,   0b111, 3);

    INIT_FLUSH_BITS(ICER_ENC_BIN_4,      0b0, 1, 1, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_4,     0b00, 2, 0, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_4,      0b1, 1, 0, 1);

    INIT_CODING_SCHEME(ICER_ENC_BIN_5,    0b00, 2,     0b1, 1);
    INIT_CODING_SCHEME(ICER_ENC_BIN_5,   0b010, 3,   0b000, 3);
    INIT_CODING_SCHEME(ICER_ENC_BIN_5,   0b110, 3,  0b0101, 4);
    INIT_CODING_SCHEME(ICER_ENC_BIN_5,   0b101, 3,  0b0100, 4);
    INIT_CODING_SCHEME(ICER_ENC_BIN_5,  0b1001, 4,  0b0111, 4);
    INIT_CODING_SCHEME(ICER_ENC_BIN_5, 0b00001, 5,  0b0010, 4);
    INIT_CODING_SCHEME(ICER_ENC_BIN_5, 0b10001, 5, 0b01100, 5);
    INIT_CODING_SCHEME(ICER_ENC_BIN_5,   0b011, 3,  0b0011, 4);
    INIT_CODING_SCHEME(ICER_ENC_BIN_5,   0b111, 3, 0b01101, 5);

    INIT_FLUSH_BITS(ICER_ENC_BIN_5,      0b0, 1, 0, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_5,     0b10, 2, 0, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_5,     0b01, 2, 1, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_5,    0b001, 3, 1, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_5,   0b0001, 4, 0, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_5,     0b11, 2, 0, 1);

    INIT_CODING_SCHEME(ICER_ENC_BIN_6,     0b1, 1,    0b01, 2);
    INIT_CODING_SCHEME(ICER_ENC_BIN_6,   0b010, 3,   0b110, 3);
    INIT_CODING_SCHEME(ICER_ENC_BIN_6,   0b110, 3,  0b1111, 4);
    INIT_CODING_SCHEME(ICER_ENC_BIN_6,   0b100, 3,   0b101, 3);
    INIT_CODING_SCHEME(ICER_ENC_BIN_6,  0b1000, 4,   0b100, 3);
    INIT_CODING_SCHEME(ICER_ENC_BIN_6, 0b10000, 5,  0b1110, 4);
    INIT_CODING_SCHEME(ICER_ENC_BIN_6, 0b00000, 5,    0b00, 2);

    INIT_FLUSH_BITS(ICER_ENC_BIN_6,     0b01, 2, 0, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_6,     0b00, 2, 1, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_6,    0b000, 3, 1, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_6,   0b0000, 4, 0, 1);

    INIT_CODING_SCHEME(ICER_ENC_BIN_7,   0b000, 3,     0b0, 1);
    INIT_CODING_SCHEME(ICER_ENC_BIN_7,   0b100, 3,   0b100, 3);
    INIT_CODING_SCHEME(ICER_ENC_BIN_7,   0b010, 3,   0b101, 3);
    INIT_CODING_SCHEME(ICER_ENC_BIN_7,   0b110, 3, 0b11110, 5);
    INIT_CODING_SCHEME(ICER_ENC_BIN_7,    0b11, 2,  0b1110, 4);
    INIT_CODING_SCHEME(ICER_ENC_BIN_7,   0b001, 3,   0b110, 3);
    INIT_CODING_SCHEME(ICER_ENC_BIN_7,   0b101, 3, 0b11111, 5);

    INIT_FLUSH_BITS(ICER_ENC_BIN_7,      0b0, 1,00, 2);
    INIT_FLUSH_BITS(ICER_ENC_BIN_7,     0b00, 2, 0, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_7,     0b10, 2, 0, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_7,      0b1, 1, 1, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_7,     0b01, 2, 0, 1);

    INIT_CODING_SCHEME(ICER_ENC_BIN_8,    0b10, 2,   0b101, 3);
    INIT_CODING_SCHEME(ICER_ENC_BIN_8,   0b100, 3,   0b100, 3);
    INIT_CODING_SCHEME(ICER_ENC_BIN_8,  0b0000, 4,     0b0, 1);
    INIT_CODING_SCHEME(ICER_ENC_BIN_8, 0b01000, 5,  0b1110, 4);
    INIT_CODING_SCHEME(ICER_ENC_BIN_8, 0b11000, 5, 0b11110, 5);
    INIT_CODING_SCHEME(ICER_ENC_BIN_8,    0b01, 2,   0b110, 3);
    INIT_CODING_SCHEME(ICER_ENC_BIN_8,    0b11, 2, 0b11111, 5);

    INIT_FLUSH_BITS(ICER_ENC_BIN_8,      0b0, 1, 1, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_8,     0b00, 2, 1, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_8,    0b000, 3, 0, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_8,   0b1000, 4, 0, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_8,      0b1, 1, 0, 1);

    return ICER_RESULT_OK;
}

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

uint16_t encode_circ_buf[2048];

int compress_partition_uint8(uint8_t *data, partition_param_typdef *params, size_t rowstride, packet_context *pkt_context,
                             output_data_buf_typedef *output_data) {
    size_t segment_w, segment_h;
    uint8_t *segment_start;
    uint16_t segment_num = 0;

    size_t partition_col_ind;
    size_t partition_row_ind = 0;

    icer_context_model_typedef context_model;
    encoder_context_typedef context;
    image_segment_typedef *seg;

    uint32_t data_in_bytes;
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

            init_context_model_vals(&context_model, pkt_context->subband_type);
            printf("t segment no: %d\n", segment_num);
            if (icer_allocate_data_packet(&seg, output_data, segment_num, pkt_context) == ICER_BYTE_QUOTA_EXCEEDED) {
                return ICER_BYTE_QUOTA_EXCEEDED;
            }
            printf("max seg out: %u\n", seg->data_length);
            init_entropy_coder_context(&context, encode_circ_buf, 2048, seg->data, seg->data_length);
            if (compress_bitplane_uint8(segment_start, segment_w, segment_h, rowstride, &context_model, &context, pkt_context) == ICER_BYTE_QUOTA_EXCEEDED) {
                data_in_bytes = context.output_ind + (context.output_bit_offset > 0);
                output_data->size_used += data_in_bytes;
                printf("exceeded seg out: %zu\n", output_data->size_used);
                return ICER_BYTE_QUOTA_EXCEEDED;
            }
            data_in_bytes = context.output_ind + (context.output_bit_offset > 0);
            seg->data_length = context.output_ind * 8 + context.output_bit_offset;
            seg->crc32 = icer_calculate_packet_crc32(seg);
            output_data->size_used += data_in_bytes;
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

            init_context_model_vals(&context_model, pkt_context->subband_type);
            printf("b segment no: %d\n", segment_num);
            if (icer_allocate_data_packet(&seg, output_data, segment_num, pkt_context) == ICER_BYTE_QUOTA_EXCEEDED) {
                return ICER_BYTE_QUOTA_EXCEEDED;
            }
            init_entropy_coder_context(&context, encode_circ_buf, 2048, seg->data, seg->data_length);
            if (compress_bitplane_uint8(segment_start, segment_w, segment_h, rowstride, &context_model, &context, pkt_context) == ICER_BYTE_QUOTA_EXCEEDED) return ICER_BYTE_QUOTA_EXCEEDED;
            data_in_bytes = context.output_ind + (context.output_bit_offset > 0);
            seg->data_length = context.output_ind * 8 + context.output_bit_offset;
            seg->crc32 = icer_calculate_packet_crc32(seg);
            output_data->size_used += data_in_bytes;
        }
        partition_row_ind += segment_h;
    }

    return ICER_RESULT_OK;
}

int compress_bitplane_uint8(uint8_t *data, size_t plane_w, size_t plane_h, size_t rowstride,
                            icer_context_model_typedef *context_model,
                            encoder_context_typedef *encoder_context,
                            packet_context* pkt_context) {
    uint8_t *pos;
    uint8_t *rowstart = data;
    int category;
    bool bit;
    uint8_t lsb = pkt_context->lsb;
    uint8_t mask = 0b1 << (lsb-1);

    printf("max out: %zu\n", encoder_context->max_output_length);

    uint8_t *h0, *h1, *v0, *v1, *d0, *d1, *d2, *d3;
    uint8_t h, v, d, tmp;
    int8_t sh0, sh1, sv0, sv1;
    uint8_t sh, sv;
    uint8_t pred_sign, actual_sign, agreement_bit;
    enum icer_pixel_contexts sign_context;
    size_t vert_bound = plane_h - 1;
    size_t hor_bound = plane_w - 1;

    uint8_t prev_plane = lsb+1;
    if (prev_plane >= 8) return ICER_BITPLANE_OUT_OF_RANGE;

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
                icer_encode_bit(encoder_context, bit, 1, 2);
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

                if (icer_encode_bit(encoder_context, bit, context_model->zero_count[cntxt], context_model->total_count[cntxt]) == ICER_BYTE_QUOTA_EXCEEDED) {
                    return ICER_BYTE_QUOTA_EXCEEDED;
                }

                context_model->total_count[cntxt]++;
                context_model->zero_count[cntxt] += (uint32_t)(!bit);
                if (context_model->total_count[cntxt] >= ICER_CONTEXT_RESCALING_CAP) {
                    context_model->total_count[cntxt] >>= 1;
                    if (context_model->zero_count[cntxt] > context_model->total_count[cntxt]) context_model->zero_count[cntxt] >>= 1;
                    else icer_ceil_div_uint32(context_model->zero_count[cntxt], 2);
                }

                if (category == ICER_CATEGORY_0 && bit) {
                    /* bit is the first magnitude bit to be encoded, thus we will encode the sign bit */
                    /* predict sign bit */
                    sh0 = 0; sh1 = 0;
                    sv0 = 0; sv1 = 0;

                    /* consider the horizontally adjacent pixels */
                    if (col > 0) sh0 = get_sign_uint8(h0, lsb);
                    if (col < hor_bound) sh1 = get_sign_uint8(h1, prev_plane);

                    /* consider the vertically adjacent pixels */
                    if (row > 0) sv0 = get_sign_uint8(v0, lsb);
                    if (row < vert_bound) sv1 = get_sign_uint8(v1, prev_plane);

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

                    if (icer_encode_bit(encoder_context, agreement_bit, context_model->zero_count[sign_context], context_model->total_count[sign_context]) == ICER_BYTE_QUOTA_EXCEEDED) {
                        return ICER_BYTE_QUOTA_EXCEEDED;
                    }

                    context_model->total_count[sign_context]++;
                    context_model->zero_count[sign_context] += (uint32_t)(agreement_bit == 0);
                    if (context_model->total_count[sign_context] >= ICER_CONTEXT_RESCALING_CAP) {
                        context_model->total_count[sign_context] >>= 1;
                        if (context_model->zero_count[sign_context] > context_model->total_count[sign_context]) context_model->zero_count[sign_context] >>= 1;
                        else icer_ceil_div_uint32(context_model->zero_count[sign_context], 2);
                    }
                }
            }

            pos++;
            /* increment position pointers of 8 adjacent pixels */
            h0++; h1++; v0++; v1++; d0++; d1++; d2++; d2++;
        }
        rowstart += rowstride;
    }
    if (flush_encode(encoder_context) == ICER_BYTE_QUOTA_EXCEEDED) {
        printf("end flush exceeded output size: %zu bytes\n", encoder_context->output_ind);
        return ICER_BYTE_QUOTA_EXCEEDED;
    }
    return ICER_RESULT_OK;
}

int decompress_bitplane_uint8(uint8_t *data, size_t plane_w, size_t plane_h, size_t rowstride,
                            icer_context_model_typedef *context_model,
                            encoder_context_typedef *encoder_context,
                            packet_context* pkt_context) {
    uint8_t *pos;
    uint8_t *rowstart = data;
    int category;
    uint8_t bit;
    uint8_t lsb = pkt_context->lsb;
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
    if (prev_plane >= 8) return ICER_BITPLANE_OUT_OF_RANGE;

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

            if (category == ICER_CATEGORY_3) {
                /* pass to uncoded bin */
                icer_decode_bit(encoder_context, &bit, 1, 2);
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

                icer_decode_bit(encoder_context, &bit, context_model->zero_count[cntxt], context_model->total_count[cntxt]);

                context_model->total_count[cntxt]++;
                context_model->zero_count[cntxt] += (uint32_t)(!bit);
                if (context_model->total_count[cntxt] >= ICER_CONTEXT_RESCALING_CAP) {
                    context_model->total_count[cntxt] >>= 1;
                    if (context_model->zero_count[cntxt] > context_model->total_count[cntxt]) context_model->zero_count[cntxt] >>= 1;
                    else icer_ceil_div_uint32(context_model->zero_count[cntxt], 2);
                }

                if (category == ICER_CATEGORY_0 && bit) {
                    /* bit is the first magnitude bit to be encoded, thus we will encode the sign bit */
                    /* predict sign bit */
                    sh0 = 0; sh1 = 0;
                    sv0 = 0; sv1 = 0;

                    /* consider the horizontally adjacent pixels */
                    if (col > 0) sh0 = get_sign_uint8(h0, lsb);
                    if (col < hor_bound) sh1 = get_sign_uint8(h1, prev_plane);

                    /* consider the vertically adjacent pixels */
                    if (row > 0) sv0 = get_sign_uint8(v0, lsb);
                    if (row < vert_bound) sv1 = get_sign_uint8(v1, prev_plane);

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

                    //todo: sign bit

                    context_model->total_count[sign_context]++;
                    context_model->zero_count[sign_context] += (uint32_t)(agreement_bit == 0);
                    if (context_model->total_count[sign_context] >= ICER_CONTEXT_RESCALING_CAP) {
                        context_model->total_count[sign_context] >>= 1;
                        if (context_model->zero_count[sign_context] > context_model->total_count[sign_context]) context_model->zero_count[sign_context] >>= 1;
                        else icer_ceil_div_uint32(context_model->zero_count[sign_context], 2);
                    }
                }
            }

            pos++;
            /* increment position pointers of 8 adjacent pixels */
            h0++; h1++; v0++; v1++; d0++; d1++; d2++; d2++;
        }
        rowstart += rowstride;
    }
    if (flush_encode(encoder_context) == ICER_BYTE_QUOTA_EXCEEDED) {
        return ICER_BYTE_QUOTA_EXCEEDED;
    }
    return ICER_RESULT_OK;
}

int icer_popbuf_while_avail(encoder_context_typedef *encoder_context) {
    uint16_t out, bits;
    //int16_t mask;
    uint16_t d, r;
    int bits_to_encode;
    while (encoder_context->used > 0 && (encoder_context->encode_buffer[encoder_context->head] & ICER_ENC_BUF_DONE_MASK)) {
        out = pop_buf(encoder_context);
        bits = out >> ICER_ENC_BUF_BITS_OFFSET;
        while (bits) {
            bits_to_encode = icer_min_int(8-encoder_context->output_bit_offset, bits);
            //mask = INT16_MIN >> (bits_to_encode-1);
            encoder_context->output_buffer[encoder_context->output_ind] |= (out & ((1 << bits_to_encode) - 1)) << encoder_context->output_bit_offset;
            out >>= bits_to_encode;
            bits -= bits_to_encode;
            r = (encoder_context->output_bit_offset + bits_to_encode) / 8;
            d = (encoder_context->output_bit_offset + bits_to_encode) % 8;
            encoder_context->output_bit_offset = d;
            if (r) {
                encoder_context->output_ind += r;
                encoder_context->output_buffer[encoder_context->output_ind] = 0;
            }
            if (encoder_context->output_ind == encoder_context->max_output_length) {
                return ICER_BYTE_QUOTA_EXCEEDED;
            }
        }
    }
    return ICER_RESULT_OK;
}

int icer_encode_bit(encoder_context_typedef *encoder_context, uint8_t bit, uint32_t zero_cnt, uint32_t total_cnt) {
    uint16_t *curr_bin;
    if (zero_cnt < (total_cnt >> 1)) {
        /*
         * we may assume that the probability of zero for each bit is contained
         * the interval [1/2, 1]
         * in the case that the probability of zero < 1/2
         * we simple invert the bit, and its associated probability
         * this is duplicated in the decoder
         */
        zero_cnt = total_cnt - zero_cnt;
        bit = bit ^ 0b1;
    }

    int bin = icer_compute_bin(zero_cnt, total_cnt);
    uint16_t prefix;
    uint16_t golomb_k;
    uint16_t bit16 = (bit != 0);

    if (encoder_context->bin_current_buf[bin] == -1) {
        encoder_context->bin_current_buf[bin] = alloc_buf(encoder_context);
        if (encoder_context->bin_current_buf[bin] == -1) {
            if (flush_encode(encoder_context) == ICER_BYTE_QUOTA_EXCEEDED) return ICER_BYTE_QUOTA_EXCEEDED;
            else encoder_context->bin_current_buf[bin] = alloc_buf(encoder_context);
        }
        curr_bin = encoder_context->encode_buffer + encoder_context->bin_current_buf[bin];
        (*curr_bin) = bin << ICER_ENC_BUF_BITS_OFFSET;
    } else curr_bin = encoder_context->encode_buffer + encoder_context->bin_current_buf[bin];

    if (bin > ICER_ENC_BIN_8) {
        /* golomb code bins */
        if (!bit16) (*curr_bin)++;
        if (bit16) {
            golomb_k = (*curr_bin) & ICER_ENC_BUF_DATA_MASK;
            (*curr_bin) = ((golomb_coders[bin].l + (golomb_k >= golomb_coders[bin].i)) << ICER_ENC_BUF_BITS_OFFSET);
            (*curr_bin) |= ((golomb_k + (golomb_k >= golomb_coders[bin].i)) & ICER_ENC_BUF_DATA_MASK);
            (*curr_bin) |= ICER_ENC_BUF_DONE_MASK;
            encoder_context->bin_current_buf[bin] = -1;
        } else if (((*curr_bin) & ICER_ENC_BUF_DATA_MASK) >= golomb_coders[bin].m) {
            (*curr_bin) = 1 << ICER_ENC_BUF_BITS_OFFSET;
            (*curr_bin) |= 1;
            (*curr_bin) |= ICER_ENC_BUF_DONE_MASK;
            encoder_context->bin_current_buf[bin] = -1;
        }
    } else if (bin != ICER_ENC_BIN_1) {
        /* custom non prefix code bins */
        (*curr_bin) |= (bit16 << encoder_context->bin_current_buf_bits[bin]);
        encoder_context->bin_current_buf_bits[bin]++;
        prefix = (*curr_bin) & ICER_ENC_BUF_DATA_MASK;
        if (custom_coding_scheme[bin][prefix].input_code_bits == encoder_context->bin_current_buf_bits[bin]) {
            (*curr_bin) = custom_coding_scheme[bin][prefix].output_code_bits << ICER_ENC_BUF_BITS_OFFSET;
            (*curr_bin) |= custom_coding_scheme[bin][prefix].output_code & ICER_ENC_BUF_DATA_MASK;
            (*curr_bin) |= ICER_ENC_BUF_DONE_MASK;
            encoder_context->bin_current_buf[bin] = -1;
            encoder_context->bin_current_buf_bits[bin] = 0;
        }

        if ((bin == ICER_ENC_BIN_5) && (prefix == 0)) {
            //printf("oops_!!!\n");
        }
    } else {
        /* uncoded bin */
        (*curr_bin) = bit16 & 0b1;
        (*curr_bin) |= (1 << ICER_ENC_BUF_BITS_OFFSET);
        (*curr_bin) |= ICER_ENC_BUF_DONE_MASK;
        encoder_context->bin_current_buf[bin] = -1;
    }

    if (icer_popbuf_while_avail(encoder_context) == ICER_BYTE_QUOTA_EXCEEDED) {
        return ICER_BYTE_QUOTA_EXCEEDED;
    }
    return ICER_RESULT_OK;
}

int icer_decode_bit(decoder_context_typedef *decoder_context, uint8_t *bit, uint32_t zero_cnt, uint32_t total_cnt) {
    uint16_t *curr_bin;
    bool inv = false;
    if (zero_cnt < (total_cnt >> 1)) {
        /*
         * we may assume that the probability of zero for each bit is contained
         * the interval [1/2, 1]
         * in the case that the probability of zero < 1/2
         * we simple invert the bit, and its associated probability
         * this is duplicated in the decoder
         */
        zero_cnt = total_cnt - zero_cnt;
        inv = true;
    }

    int bin = icer_compute_bin(zero_cnt, total_cnt);

    if (bin > ICER_ENC_BIN_8) {
        /* golomb code bins */

    } else if (bin != ICER_ENC_BIN_1) {
        /* custom non prefix code bins */

    } else {
        /* uncoded bin */

    }

    //if (icer_popbuf_while_avail(encoder_context) == ICER_BYTE_QUOTA_EXCEEDED) return ICER_BYTE_QUOTA_EXCEEDED;
    return ICER_RESULT_OK;
}

int flush_encode(encoder_context_typedef *encoder_context) {
    //printf("to flush: %3zu codewords\n", encoder_context->used);
    uint16_t *first = encoder_context->encode_buffer + encoder_context->head;
    custom_flush_typedef *flush;
    uint16_t prefix;
    if (((*first) & ICER_ENC_BUF_DONE_MASK) == 0) {
        uint8_t bin = (*first) >> ICER_ENC_BUF_BITS_OFFSET;

        uint16_t golomb_k;

        if (bin > ICER_ENC_BIN_8) {
            /* golomb code bins */
            golomb_k = (*first) & ICER_ENC_BUF_DATA_MASK;
            if (golomb_k == golomb_coders[bin].m - 1) {
                *first = 1 << ICER_ENC_BUF_BITS_OFFSET;
                *first |= 1;
                *first |= ICER_ENC_BUF_DONE_MASK;
            } else {
                *first = ((golomb_coders[bin].l + (golomb_k >= golomb_coders[bin].i)) << ICER_ENC_BUF_BITS_OFFSET);
                *first |= ((golomb_k + (golomb_k >= golomb_coders[bin].i)) & ICER_ENC_BUF_DATA_MASK);
                *first |= ICER_ENC_BUF_DONE_MASK;
            }
            encoder_context->bin_current_buf[bin] = -1;
        } else if (bin != ICER_ENC_BIN_1) {
            /* custom non prefix code bins */
            flush = &(custom_code_flush_bits[bin][(*first) & ICER_ENC_BUF_DATA_MASK][encoder_context->bin_current_buf_bits[bin]]);
            (*first) |= flush->flush_bit << encoder_context->bin_current_buf_bits[bin];
            encoder_context->bin_current_buf_bits[bin] += flush->flush_bit_numbers;

            prefix = (*first) & ICER_ENC_BUF_DATA_MASK;

            (*first) = custom_coding_scheme[bin][prefix].output_code_bits << ICER_ENC_BUF_BITS_OFFSET;
            (*first) |= custom_coding_scheme[bin][prefix].output_code & ICER_ENC_BUF_DATA_MASK;
            (*first) |= ICER_ENC_BUF_DONE_MASK;
            encoder_context->bin_current_buf[bin] = -1;
            encoder_context->bin_current_buf_bits[bin] = 0;
        } else {
            /* uncoded bin */
            // this should never happen?
            printf("oopsies its broken\n");
        }
    }

    if (icer_popbuf_while_avail(encoder_context) == ICER_BYTE_QUOTA_EXCEEDED) return ICER_BYTE_QUOTA_EXCEEDED;
    return ICER_RESULT_OK;
}

void init_context_model_vals(icer_context_model_typedef* context_model, enum icer_subband_types subband_type) {
    context_model->subband_type = subband_type;
    for (size_t i = 0;i <= ICER_CONTEXT_MAX;i++) {
        context_model->zero_count[i] = ICER_DEFAULT_CONTEXT_ZERO_COUNT;
        context_model->total_count[i] = ICER_DEFAULT_CONTEXT_TOTAL_COUNT;
    }
}

void init_entropy_coder_context(encoder_context_typedef *encoder_context, uint16_t *encode_buffer, size_t buffer_length, uint8_t *encoder_out, size_t enc_out_max) {
    encoder_context->max_output_length = enc_out_max;
    encoder_context->output_buffer = encoder_out;

    encoder_context->buffer_length = buffer_length;
    encoder_context->encode_buffer = encode_buffer;

    encoder_context->head = 0;
    encoder_context->tail = 0;
    encoder_context->used = 0;

    for (size_t it = 0;it < ICER_ENCODER_BIN_MAX+1;it++) {
        encoder_context->bin_current_buf[it] = -1;
        encoder_context->bin_current_buf_bits[it] = 0;
    }

    encoder_context->output_ind = 0;
    encoder_context->output_bit_offset = 0;
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

uint32_t icer_calculate_packet_crc32(image_segment_typedef *pkt) {
    return crc32buf((char*) pkt, 10);
}

int icer_allocate_data_packet(image_segment_typedef **pkt, output_data_buf_typedef *output_data, uint8_t segment_num, packet_context *context) {
    size_t buf_len = output_data->size_allocated - output_data->size_used;
    if (buf_len < 16) {
        return ICER_BYTE_QUOTA_EXCEEDED;
    }
    (*pkt) = (image_segment_typedef *) (output_data->data_start + output_data->size_used);
    (*pkt)->preamble = ICER_PACKET_PREAMBLE;
    (*pkt)->subband_number = context->subband_number;
    (*pkt)->subband_type = context->subband_type;
    (*pkt)->segment_number = segment_num;
    (*pkt)->lsb = context->lsb;
    (*pkt)->ll_mean_val = context->ll_mean_val;
    (*pkt)->crc32 = 0;

    output_data->size_used += 16;
    buf_len -= 16;

    // store max data length first
    (*pkt)->data_length = buf_len;
    (*pkt)->data = (uint8_t *) (output_data->data_start + output_data->size_used);

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

void icer_init_output_struct(output_data_buf_typedef *out, uint8_t *data, size_t len) {
    out->size_used = 0;
    out->data_start = data;
    out->size_allocated = len;
}

static inline uint16_t pop_buf(encoder_context_typedef *cntxt) {
    if (cntxt->used > 0) cntxt->used--;
    uint16_t res = cntxt->encode_buffer[cntxt->head];
    cntxt->head = (cntxt->head + 1) % cntxt->buffer_length;
    return res;
}

static inline int16_t alloc_buf(encoder_context_typedef *cntxt) {
    if (cntxt->used+1 >= cntxt->buffer_length) return -1;
    cntxt->used++;
    int16_t ind = (int16_t)cntxt->tail;
    cntxt->tail = (cntxt->tail+1) % cntxt->buffer_length;
    return ind;
}

int icer_compute_bin(uint32_t zero_cnt, uint32_t total_cnt) {
    uint32_t comp = zero_cnt * ICER_BIN_PROBABILITY_DENOMINATOR;
    for (int16_t bin = ICER_ENCODER_BIN_MAX;bin > ICER_ENC_BIN_1;bin--) {
        if (comp >= total_cnt * icer_bin_probability_cutoffs[bin-1]) {
            return bin;
        }
    }
    return ICER_ENC_BIN_1;
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

/* optimizes into single division, as per https://stackoverflow.com/questions/46265403/fast-floor-of-a-signed-integer-division-in-c-c */
static inline int16_t icer_floor_div_int16(int16_t a, int16_t b) {
    int16_t d = a / b;
    int16_t r = a % b;
    return r ? (d - ((a < 0) ^ (b < 0))) : d;
}

static inline size_t icer_floor_div_size_t(size_t a, size_t b) {
    return a / b;
}

static inline int16_t icer_ceil_div_int16(int16_t a, int16_t b) {
    int16_t d = a / b;
    int16_t r = a % b;
    return r ? d + 1 : d;
}

static inline uint32_t icer_ceil_div_uint32(uint32_t a, uint32_t b) {
    uint32_t d = a / b;
    uint32_t r = a % b;
    return r ? d + 1 : d;
}

static inline size_t icer_ceil_div_size_t(size_t a, size_t b) {
    size_t d = a / b;
    size_t r = a % b;
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
