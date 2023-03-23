//
// Created by linyi on 19/3/2023.
//

#include "../inc/icer.h"

int icer_compute_bin(uint32_t zero_cnt, uint32_t total_cnt) {
    uint32_t comp = zero_cnt * ICER_BIN_PROBABILITY_DENOMINATOR;
    for (int16_t bin = ICER_ENCODER_BIN_MAX;bin > ICER_ENC_BIN_1;bin--) {
        if (comp >= total_cnt * icer_bin_probability_cutoffs[bin-1]) {
            return bin;
        }
    }
    return ICER_ENC_BIN_1;
}

uint32_t icer_calculate_packet_crc32(icer_image_segment_typedef *pkt) {
    return crc32buf((char*)pkt, sizeof(icer_image_segment_typedef) - 4);
}

uint32_t icer_calculate_segment_crc32(icer_image_segment_typedef *pkt) {
    return crc32buf((char*)pkt + sizeof(icer_image_segment_typedef), icer_ceil_div_uint32(pkt->data_length, 8));
}

