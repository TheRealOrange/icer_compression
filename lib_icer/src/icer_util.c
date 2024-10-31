//
// Created by linyi on 19/3/2023.
//

#include "icer.h"

#ifndef CRC32BUF_FUNCTION
#include "crc.h"
#define CRC32BUF_FUNCTION(x, y) crc32buf(x, y)
#endif

#ifndef USER_PROVIDED_BUFFERS

#ifdef USE_UINT8_FUNCTIONS
#ifdef USE_ENCODE_FUNCTIONS
icer_packet_context icer_packets[ICER_MAX_PACKETS];
icer_image_segment_typedef *icer_rearrange_segments_8[ICER_CHANNEL_MAX + 1][ICER_MAX_DECOMP_STAGES + 1][ICER_SUBBAND_MAX + 1][7][ICER_MAX_SEGMENTS + 1];
#endif

#ifdef USE_DECODE_FUNCTIONS
icer_image_segment_typedef *icer_reconstruct_data_8[ICER_CHANNEL_MAX + 1][ICER_MAX_DECOMP_STAGES + 1][ICER_SUBBAND_MAX + 1][ICER_MAX_SEGMENTS + 1][7];
#endif
#endif

#ifdef USE_UINT16_FUNCTIONS
#ifdef USE_ENCODE_FUNCTIONS
icer_packet_context icer_packets_16[ICER_MAX_PACKETS_16];
icer_image_segment_typedef *icer_rearrange_segments_16[ICER_CHANNEL_MAX + 1][ICER_MAX_DECOMP_STAGES + 1][ICER_SUBBAND_MAX + 1][15][ICER_MAX_SEGMENTS + 1];
#endif

#ifdef USE_DECODE_FUNCTIONS
icer_image_segment_typedef *icer_reconstruct_data_16[ICER_CHANNEL_MAX + 1][ICER_MAX_DECOMP_STAGES + 1][ICER_SUBBAND_MAX + 1][ICER_MAX_SEGMENTS + 1][15];
#endif
#endif

#endif

int icer_init_output_struct(icer_output_data_buf_typedef *out, uint8_t *data, size_t buf_len, size_t byte_quota) {
    if (byte_quota * 2 > buf_len) return ICER_OUTPUT_BUF_TOO_SMALL;
    out->size_used = 0;
    out->data_start = data;
    out->size_allocated = byte_quota;
    out->rearrange_start = data + byte_quota;
    return ICER_RESULT_OK;
}

/* compute which bin of the interleaved entropy coder to place a bit to be encoded based of the probability cutoffs of each bin */
int icer_compute_bin(uint32_t zero_cnt, uint32_t total_cnt) {
    uint32_t comp = zero_cnt * ICER_BIN_PROBABILITY_DENOMINATOR;
    for (int16_t bin = ICER_ENCODER_BIN_MAX;bin > ICER_ENC_BIN_1;bin--) {
        if (comp >= total_cnt * icer_bin_probability_cutoffs[bin-1]) {
            return bin;
        }
    }
    return ICER_ENC_BIN_1;
}

/* calculates the crc32 for the header of each segment of the image */
uint32_t icer_calculate_packet_crc32(icer_image_segment_typedef *pkt) {
    return crc32buf((char*)pkt, sizeof(icer_image_segment_typedef) - 4);
}

/* calculates the crc32 for the data portion of each segment of the image */
uint32_t icer_calculate_segment_crc32(icer_image_segment_typedef *pkt) {
    return crc32buf((char*)pkt + sizeof(icer_image_segment_typedef), icer_ceil_div_uint32(pkt->data_length, 8));
}

#ifdef USE_DECODE_FUNCTIONS
#ifdef USE_UINT8_FUNCTIONS
void icer_remove_negative_uint8(uint8_t * const image, size_t image_w, size_t image_h) {
    size_t data_length = image_w * image_h;
    int8_t *data_end = (int8_t *)(image + data_length);
    for (int8_t *pixel = (int8_t *)image;pixel < data_end;pixel++) {
        if (*pixel < 0) {
            *pixel = 0;
        }
    }
}
#endif

#ifdef USE_UINT16_FUNCTIONS
void icer_remove_negative_uint16(uint16_t * const image, size_t image_w, size_t image_h) {
    size_t data_length = image_w * image_h;
    int16_t *data_end = (int16_t *)(image + data_length);
    for (int16_t *pixel = (int16_t *)image;pixel < data_end;pixel++) {
        if (*pixel < 0) {
            *pixel = 0;
        }
    }
}
#endif
#endif