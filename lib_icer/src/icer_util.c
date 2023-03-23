//
// Created by linyi on 19/3/2023.
//

#include "../inc/icer.h"

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

