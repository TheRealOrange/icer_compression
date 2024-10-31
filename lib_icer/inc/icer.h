#ifndef ICER_COMPRESSION_ICER_H
#define ICER_COMPRESSION_ICER_H

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef _MSC_VER
#include <intrin.h>

#define __builtin_popcount __popcnt
#define __builtin_clz __clz

uint32_t __inline __clz( uint32_t value ) {
    uint32_t leading_zero = 0;

    if ( _BitScanReverse( &leading_zero, value ) ) {
       return 31 - leading_zero;
    } else {
         // Same remarks as above
         return 32;
    }
}
#endif

#define ICER_CIRC_BUF_SIZE 2048
#define MAX_K 12
#ifndef ICER_MAX_SEGMENTS
#define ICER_MAX_SEGMENTS 32
#endif
#ifndef ICER_MAX_DECOMP_STAGES
#define ICER_MAX_DECOMP_STAGES 6
#endif
#ifndef ICER_MAX_PACKETS
#define ICER_MAX_PACKETS 300
#endif
#ifndef ICER_MAX_PACKETS_16
#define ICER_MAX_PACKETS_16 800
#endif
#ifndef ICER_BITPLANES_TO_COMPRESS_8
#define ICER_BITPLANES_TO_COMPRESS_8 7
#endif
#ifndef ICER_BITPLANES_TO_COMPRESS_16
#define ICER_BITPLANES_TO_COMPRESS_16 9
#endif

//#define USER_PROVIDED_BUFFERS
/*
 * if the user decides to specify explicitly where to place the buffers used during encoding and
 * decoding, the user should define the above line and allocate memory for the buffers listed below
 *
 * for encoding uint8
 * icer_packet_context icer_packets[ICER_MAX_PACKETS];
 * icer_image_segment_typedef *icer_rearrange_segments_8[ICER_CHANNEL_MAX + 1][ICER_MAX_DECOMP_STAGES + 1][ICER_SUBBAND_MAX + 1][7][ICER_MAX_SEGMENTS + 1];
 *
 * for decoding uint8
 * icer_image_segment_typedef *icer_reconstruct_data_8[ICER_CHANNEL_MAX + 1][ICER_MAX_DECOMP_STAGES + 1][ICER_SUBBAND_MAX + 1][ICER_MAX_SEGMENTS + 1][7];
 *
 * for encoding uint16
 * icer_packet_context icer_packets_16[ICER_MAX_PACKETS_16];
 * icer_image_segment_typedef *icer_rearrange_segments_16[ICER_CHANNEL_MAX + 1][ICER_MAX_DECOMP_STAGES + 1][ICER_SUBBAND_MAX + 1][15][ICER_MAX_SEGMENTS + 1];
 *
 * for decoding uint16
 * icer_image_segment_typedef *icer_reconstruct_data_16[ICER_CHANNEL_MAX + 1][ICER_MAX_DECOMP_STAGES + 1][ICER_SUBBAND_MAX + 1][ICER_MAX_SEGMENTS + 1][15];
 *
 * common to all encoding
 * uint16_t icer_encode_circ_buf[ICER_CIRC_BUF_SIZE];
 */

//#define CRC32BUF_FUNCTION(x, y) xxx
/*
 *
 * if the user decides to use a custom provided crc32 implementation (perhaps in hardware) the
 * user should define the above line with the function as specified below
 *
 * uint32_t crc32buf(char *buf, size_t len)
 *
 * override function with a custom defined function
 */

#if !defined(USE_UINT8_FUNCTIONS) && !defined(USE_UINT16_FUNCTIONS)
#define USE_UINT8_FUNCTIONS
#define USE_UINT16_FUNCTIONS
#endif

#if !defined(USE_DECODE_FUNCTIONS) && !defined(USE_ENCODE_FUNCTIONS)
#define USE_DECODE_FUNCTIONS
#define USE_ENCODE_FUNCTIONS
#endif

enum icer_status {
    ICER_RESULT_OK = 0,
    ICER_INTEGER_OVERFLOW = -1,
    ICER_OUTPUT_BUF_TOO_SMALL = -2,
    ICER_TOO_MANY_SEGMENTS = -3,
    ICER_TOO_MANY_STAGES = -4,
    ICER_BYTE_QUOTA_EXCEEDED = -5,
    ICER_BITPLANE_OUT_OF_RANGE = -6,
    ICER_DECODER_OUT_OF_DATA = -7,
    ICER_DECODED_INVALID_DATA = -8,
    ICER_PACKET_COUNT_EXCEEDED = -9,
    ICER_FATAL_ERROR = -10,
    ICER_INVALID_INPUT = -11,
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

extern const int16_t icer_wavelet_filter_parameters[][4];

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

extern const uint8_t icer_context_table_ll_lh_hl[3][3][5];

extern const uint8_t icer_context_table_hh[5][5];

extern const uint8_t icer_sign_context_table[5][5];

/* 1 is negative, 0 is positive */
extern const uint8_t icer_sign_prediction_table[5][5];

enum icer_pixel_categories {
    ICER_CATEGORY_0 = 0,
    ICER_CATEGORY_1,
    ICER_CATEGORY_2,
    ICER_CATEGORY_3 = 3
};

#define ICER_SUBBAND_MAX 3
enum icer_subband_types {
    ICER_SUBBAND_LL = 0,
    ICER_SUBBAND_HL,
    ICER_SUBBAND_LH,
    ICER_SUBBAND_HH = ICER_SUBBAND_MAX
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

extern const uint32_t icer_bin_probability_cutoffs[ICER_ENCODER_BIN_MAX+1];

#define ICER_BIN_PROBABILITY_DENOMINATOR 65536

extern const int32_t icer_bin_coding_scheme[ICER_ENCODER_BIN_MAX+1];

typedef struct {
    uint8_t input_code_bits;
    uint8_t output_code_bits;
    uint8_t output_code;
} icer_custom_code_typedef;

typedef struct {
    uint8_t flush_bit;
    uint8_t flush_bit_numbers;
} icer_custom_flush_typedef;

#define CUSTOM_CODING_MAX_LOOKUP 32
#ifdef USE_ENCODE_FUNCTIONS
extern icer_custom_code_typedef icer_custom_coding_scheme[ICER_ENCODER_BIN_MAX + 1][CUSTOM_CODING_MAX_LOOKUP];
#endif

#ifdef USE_DECODE_FUNCTIONS
extern icer_custom_code_typedef icer_custom_decode_scheme[ICER_ENCODER_BIN_MAX + 1][CUSTOM_CODING_MAX_LOOKUP];
#endif

#define CUSTOM_CODE_FLUSH_MAX_LOOKUP 8
#define MAX_NUM_BITS_BEFORE_FLUSH    5
extern icer_custom_flush_typedef icer_custom_code_flush_bits[ICER_ENCODER_BIN_MAX + 1][CUSTOM_CODE_FLUSH_MAX_LOOKUP + 1][MAX_NUM_BITS_BEFORE_FLUSH + 1];

typedef struct {
    uint16_t m;
    uint16_t l;
    uint16_t i;
} icer_golomb_code_typedef;

extern icer_golomb_code_typedef icer_golomb_coders[ICER_ENCODER_BIN_MAX + 1];

typedef struct {
    uint8_t decomp_level;
    uint8_t subband_type;
    uint8_t ll_mean_val;
    uint8_t lsb;
    uint64_t priority;
    size_t image_w;
    size_t image_h;
    uint8_t channel;
} icer_packet_context;

#ifdef USE_UINT8_FUNCTIONS
extern icer_packet_context icer_packets[ICER_MAX_PACKETS];
#endif

#ifdef USE_UINT16_FUNCTIONS
extern icer_packet_context icer_packets_16[ICER_MAX_PACKETS_16];
#endif

#define ICER_PACKET_PREAMBLE 0x605B
#define ICER_SEGMENT_LSB_MASK 0x0f
#define ICER_SEGMENT_CHANNEL_MASK 0xf0
#define ICER_GET_LSB_MACRO(x) ((x) & ICER_SEGMENT_LSB_MASK)
#define ICER_GET_CHANNEL_MACRO(x) (((x) & ICER_SEGMENT_CHANNEL_MASK) >> 4)
#define ICER_SET_CHANNEL_MACRO(x) ((x) << 4)

typedef struct {
    uint16_t preamble;
    uint16_t ll_mean_val;
    uint8_t decomp_level;
    uint8_t subband_type;
    uint8_t segment_number;
    uint8_t lsb_chan;
    uint32_t image_w;
    uint32_t image_h;
    uint32_t data_length; // store data length in bits for the decoder
    uint32_t data_crc32;
    uint32_t crc32;
} icer_image_segment_typedef;

typedef struct {
    size_t size_used;
    size_t size_allocated;
    uint8_t *data_start;
    uint8_t *rearrange_start;
} icer_output_data_buf_typedef;

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
} icer_encoder_context_typedef;

#define ICER_DECODER_BIT_BIN_MAX 30

typedef struct {
    size_t decoded_words;
    size_t encode_ind;
    uint8_t encode_bit_offset;
    uint32_t encoded_bits_total;
    uint32_t decoded_bits_total;
    uint8_t *encoded_words;
    uint32_t bin_buf[ICER_ENCODER_BIN_MAX+1][ICER_DECODER_BIT_BIN_MAX];
    int32_t bin_bits[ICER_ENCODER_BIN_MAX+1];
    size_t bin_decode_index[ICER_ENCODER_BIN_MAX+1];
} icer_decoder_context_typedef;

#define ICER_CHANNEL_MIN 0
#define ICER_CHANNEL_MAX 2
enum icer_color_channels {
    ICER_CHANNEL_Y = ICER_CHANNEL_MIN,
    ICER_CHANNEL_U = 1,
    ICER_CHANNEL_V = ICER_CHANNEL_MAX
};

#ifdef USE_UINT8_FUNCTIONS
#ifdef USE_DECODE_FUNCTIONS
extern icer_image_segment_typedef *icer_reconstruct_data_8[ICER_CHANNEL_MAX + 1][ICER_MAX_DECOMP_STAGES + 1][ICER_SUBBAND_MAX + 1][ICER_MAX_SEGMENTS + 1][7];
#endif

#ifdef USE_ENCODE_FUNCTIONS
extern icer_image_segment_typedef *icer_rearrange_segments_8[ICER_CHANNEL_MAX + 1][ICER_MAX_DECOMP_STAGES + 1][ICER_SUBBAND_MAX + 1][7][ICER_MAX_SEGMENTS + 1];
#endif
#endif

#ifdef USE_UINT16_FUNCTIONS
#ifdef USE_DECODE_FUNCTIONS
extern icer_image_segment_typedef *icer_reconstruct_data_16[ICER_CHANNEL_MAX + 1][ICER_MAX_DECOMP_STAGES + 1][ICER_SUBBAND_MAX + 1][ICER_MAX_SEGMENTS + 1][15];
#endif

#ifdef USE_ENCODE_FUNCTIONS
extern icer_image_segment_typedef *icer_rearrange_segments_16[ICER_CHANNEL_MAX + 1][ICER_MAX_DECOMP_STAGES + 1][ICER_SUBBAND_MAX + 1][15][ICER_MAX_SEGMENTS + 1];
#endif
#endif

int icer_init(void);

#ifdef USE_DECODE_FUNCTIONS
void icer_init_decodescheme(void);
int icer_get_image_dimensions(const uint8_t *datastream, size_t data_length, size_t *image_w, size_t *image_h);
#endif

#ifdef USE_ENCODE_FUNCTIONS
void icer_init_codingscheme(void);
#endif

void icer_init_flushbits(void);
void icer_init_golombcoder(void);

#ifdef USE_UINT8_FUNCTIONS
#ifdef USE_ENCODE_FUNCTIONS
int icer_compress_image_uint8(uint8_t *image, size_t image_w, size_t image_h, uint8_t stages, enum icer_filter_types filt,
                              uint8_t segments, icer_output_data_buf_typedef *output_data);
int icer_compress_image_yuv_uint8(uint8_t *y_channel, uint8_t *u_channel, uint8_t *v_channel, size_t image_w,
                                  size_t image_h, uint8_t stages, enum icer_filter_types filt,
                                  uint8_t segments, icer_output_data_buf_typedef *const output_data);

int icer_wavelet_transform_stages_uint8(uint8_t *image, size_t image_w, size_t image_h, uint8_t stages, enum icer_filter_types filt);

int icer_wavelet_transform_2d_uint8(uint8_t *image, size_t image_w, size_t image_h, size_t rowstride, enum icer_filter_types filt);
int icer_wavelet_transform_1d_uint8(uint8_t *data, size_t N, size_t stride, enum icer_filter_types filt);

int icer_compress_partition_uint8(const uint8_t *data, partition_param_typdef *params, size_t rowstride,
                                  icer_packet_context *pkt_context, icer_output_data_buf_typedef *output_data,
                                  icer_image_segment_typedef *segments_encoded[]);
int icer_compress_bitplane_uint8(const uint8_t *data, size_t plane_w, size_t plane_h, size_t rowstride,
                                 icer_context_model_typedef *context_model,
                                 icer_encoder_context_typedef *encoder_context,
                                 const icer_packet_context *pkt_context);
#endif

#ifdef USE_DECODE_FUNCTIONS
int icer_decompress_image_uint8(uint8_t *image, size_t *image_w, size_t *image_h, size_t image_bufsize, const uint8_t *datastream,
                                size_t data_length, uint8_t stages, enum icer_filter_types filt, uint8_t segments);
int icer_decompress_image_yuv_uint8(uint8_t *y_channel, uint8_t *u_channel, uint8_t *v_channel, size_t *image_w,
                                    size_t *image_h, size_t image_bufsize, const uint8_t *datastream,
                                    size_t data_length, uint8_t stages, enum icer_filter_types filt,
                                    uint8_t segments);

int icer_inverse_wavelet_transform_stages_uint8(uint8_t *image, size_t image_w, size_t image_h, uint8_t stages, enum icer_filter_types filt);

int icer_inverse_wavelet_transform_2d_uint8(uint8_t *image, size_t image_w, size_t image_h, size_t rowstride, enum icer_filter_types filt);
int icer_inverse_wavelet_transform_1d_uint8(uint8_t *data, size_t N, size_t stride, enum icer_filter_types filt);

int icer_decompress_partition_uint8(uint8_t *const data, const partition_param_typdef *params, size_t rowstride,
                                    icer_image_segment_typedef *seg[][7]);
int icer_decompress_bitplane_uint8(uint8_t *const data, size_t plane_w, size_t plane_h, size_t rowstride,
                                   icer_context_model_typedef *context_model,
                                   icer_decoder_context_typedef *decoder_context,
                                   const icer_packet_context *pkt_context);
void icer_remove_negative_uint8(uint8_t * const image, size_t image_w, size_t image_h);
#endif

void icer_reverse_uint8(uint8_t *data, size_t start, size_t end, size_t stride);

void icer_interleave_uint8(uint8_t *data, size_t len, size_t stride);
void icer_deinterleave_uint8(uint8_t *data, size_t len, size_t stride);


void icer_to_sign_magnitude_int8(uint8_t *data, size_t len);
void icer_from_sign_magnitude_int8(uint8_t *data, size_t len);
#endif

#ifdef USE_UINT16_FUNCTIONS
#ifdef USE_ENCODE_FUNCTIONS
int icer_compress_image_uint16(uint16_t *image, size_t image_w, size_t image_h, uint8_t stages, enum icer_filter_types filt,
                              uint8_t segments, icer_output_data_buf_typedef *output_data);
int icer_compress_image_yuv_uint16(uint16_t *y_channel, uint16_t *u_channel, uint16_t *v_channel, size_t image_w,
                                   size_t image_h, uint8_t stages, enum icer_filter_types filt,
                                   uint8_t segments, icer_output_data_buf_typedef *output_data);

int icer_wavelet_transform_stages_uint16(uint16_t *image, size_t image_w, size_t image_h, uint8_t stages, enum icer_filter_types filt);

int icer_wavelet_transform_2d_uint16(uint16_t *image, size_t image_w, size_t image_h, size_t rowstride, enum icer_filter_types filt);
int icer_wavelet_transform_1d_uint16(uint16_t *data, size_t N, size_t stride, enum icer_filter_types filt);

int icer_compress_partition_uint16(const uint16_t *data, const partition_param_typdef *params, size_t rowstride,
                                   const icer_packet_context *pkt_context, icer_output_data_buf_typedef *output_data,
                                   icer_image_segment_typedef *segments_encoded[]);
int icer_compress_bitplane_uint16(const uint16_t *data, size_t plane_w, size_t plane_h, size_t rowstride,
                                  icer_context_model_typedef *context_model,
                                  icer_encoder_context_typedef *encoder_context,
                                  const icer_packet_context *pkt_context);
#endif

#ifdef USE_DECODE_FUNCTIONS
int icer_decompress_image_uint16(uint16_t *image, size_t *image_w, size_t *image_h, size_t image_bufsize, const uint8_t *datastream,
                                 size_t data_length, uint8_t stages, enum icer_filter_types filt, uint8_t segments);
int icer_decompress_image_yuv_uint16(uint16_t * y_channel, uint16_t * u_channel, uint16_t * v_channel, size_t *image_w,
                                     size_t *image_h, size_t image_bufsize, const uint8_t *datastream,
                                     size_t data_length, uint8_t stages, enum icer_filter_types filt,
                                     uint8_t segments);

int icer_inverse_wavelet_transform_stages_uint16(uint16_t *image, size_t image_w, size_t image_h, uint8_t stages, enum icer_filter_types filt);

int icer_inverse_wavelet_transform_2d_uint16(uint16_t *image, size_t image_w, size_t image_h, size_t rowstride, enum icer_filter_types filt);
int icer_inverse_wavelet_transform_1d_uint16(uint16_t *data, size_t N, size_t stride, enum icer_filter_types filt);

int icer_decompress_partition_uint16(uint16_t *data, const partition_param_typdef *params, size_t rowstride,
                                     icer_image_segment_typedef *seg[][15]);
int icer_decompress_bitplane_uint16(uint16_t *const data, size_t plane_w, size_t plane_h, size_t rowstride,
                                    icer_context_model_typedef *context_model,
                                    icer_decoder_context_typedef *decoder_context,
                                    const icer_packet_context *pkt_context);
void icer_remove_negative_uint16(uint16_t * const image, size_t image_w, size_t image_h);
#endif

void icer_reverse_uint16(uint16_t *data, size_t start, size_t end, size_t stride);

void icer_interleave_uint16(uint16_t *data, size_t len, size_t stride);
void icer_deinterleave_uint16(uint16_t *data, size_t len, size_t stride);

void icer_to_sign_magnitude_int16(uint16_t *data, size_t len);
void icer_from_sign_magnitude_int16(uint16_t *data, size_t len);
#endif

#ifdef USE_ENCODE_FUNCTIONS
/* inside encoding.c */
extern uint16_t icer_encode_circ_buf[ICER_CIRC_BUF_SIZE];

void icer_init_entropy_coder_context(icer_encoder_context_typedef *encoder_context, uint16_t *encode_buffer, size_t buffer_length, uint8_t *encoder_out, size_t enc_out_max);
int icer_encode_bit(icer_encoder_context_typedef *encoder_context, uint8_t bit, uint32_t zero_cnt, uint32_t total_cnt);
int icer_popbuf_while_avail(icer_encoder_context_typedef *encoder_context);
int icer_flush_encode(icer_encoder_context_typedef *encoder_context);
int icer_allocate_data_packet(icer_image_segment_typedef **pkt, icer_output_data_buf_typedef * const output_data, uint8_t segment_num, const icer_packet_context *context);
#endif

#ifdef USE_DECODE_FUNCTIONS
/* inside decoding.c */
void icer_init_entropy_decoder_context(icer_decoder_context_typedef *decoder_context, uint8_t *encoded_words, size_t encoded_bits);
void icer_push_bin_bits(icer_decoder_context_typedef *decoder_context, uint8_t bin, uint16_t bits, uint16_t num_bits);
int icer_get_bit_from_codeword(icer_decoder_context_typedef *decoder_context, uint8_t bits);
int icer_get_bits_from_codeword(icer_decoder_context_typedef *decoder_context, uint8_t bits);
int icer_pop_bits_from_codeword(icer_decoder_context_typedef *decoder_context, uint8_t bits);
int icer_decode_bit(icer_decoder_context_typedef *decoder_context, uint8_t *bit, uint32_t zero_cnt, uint32_t total_cnt);
#endif

int icer_find_packet_in_bytestream(icer_image_segment_typedef **seg, const uint8_t *datastream, size_t data_length, size_t * offset);

uint8_t icer_find_k(size_t len);
size_t icer_get_dim_n_low_stages(size_t dim, uint8_t stages);
size_t icer_get_dim_n_high_stages(size_t dim, uint8_t stages);

void icer_init_context_model_vals(icer_context_model_typedef* context_model, enum icer_subband_types subband_type);
int icer_generate_partition_parameters(partition_param_typdef *params, size_t ll_w, size_t ll_h, uint16_t segments);

uint32_t icer_calculate_packet_crc32(icer_image_segment_typedef *pkt);
uint32_t icer_calculate_segment_crc32(icer_image_segment_typedef *pkt);

int icer_init_output_struct(icer_output_data_buf_typedef *out, uint8_t *data, size_t buf_len, size_t byte_quota);
int icer_compute_bin(uint32_t zero_cnt, uint32_t total_cnt);

static inline unsigned icer_pow_uint(unsigned base, unsigned exp);

static inline int16_t icer_floor_div_int16(int16_t a, int16_t b);
static inline int32_t icer_floor_div_int32(int32_t a, int32_t b);
static inline size_t icer_floor_div_size_t(size_t a, size_t b);
static inline int16_t icer_ceil_div_int16(int16_t a, int16_t b);
static inline uint32_t icer_ceil_div_uint32(uint32_t a, uint32_t b);
static inline size_t icer_ceil_div_size_t(size_t a, size_t b);

static inline int icer_max_int(int a, int b);
static inline unsigned icer_max_uint(unsigned a, unsigned b);

static inline int icer_min_int(int a, int b);

static inline void icer_reverse_bits(uint16_t *bits, uint8_t num);

extern size_t icer_slice_lengths[MAX_K];

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

static inline int32_t icer_floor_div_int32(int32_t a, int32_t b) {
    int32_t d = a / b;
    int32_t r = a % b;
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

static inline void icer_reverse_bits(uint16_t *bits, uint8_t num) {
    uint16_t reversed = 0;
    for (int b = 0; b < num;b++) {
        reversed <<= 1;
        reversed |= (*bits) & 1;
        (*bits) >>= 1;
    }
    (*bits) = reversed;
}


#endif //ICER_COMPRESSION_ICER_H
