//
// Created by linyi on 19/3/2023.
//

#include "icer.h"

int icer_get_bits_from_codeword(decoder_context_typedef *decoder_context, uint8_t bits) {
    int num = 0;
    int bits_to_decode= 0;
    uint8_t bitoffset = decoder_context->encode_bit_offset;
    size_t ind = decoder_context->encode_ind;
    while (bits) {
        bits_to_decode = icer_min_int(bitoffset, bits);

        if (decoder_context->decoded_bits_total + bits_to_decode > decoder_context->encoded_bits_total) {
            return ICER_DECODER_OUT_OF_DATA;
        }

        bits -= bits_to_decode;
        num |= (decoder_context->encoded_words[ind] & ((1 << bits_to_decode) - 1)) << bits;
        bitoffset -= bits_to_decode;
        if (bitoffset == 0) {
            ind++;
            bitoffset = 7;
        }
    }
    return num;
}

int icer_decode_bit(decoder_context_typedef *decoder_context, uint8_t *bit, uint32_t zero_cnt, uint32_t total_cnt) {
    bool inv = false, b;
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

    if (decoder_context->bin_bits[bin] <= 0) {
        /* ran out of bits in the bit, time to process a new codeword */
        if (bin > ICER_ENC_BIN_8) {
            /* golomb code bins */

        } else if (bin != ICER_ENC_BIN_1) {
            /* custom non prefix code bins */

        } else {
            /* uncoded bin */

        }
    }

    b = (decoder_context->bin_buf[bin][decoder_context->bin_decode_index[bin]] & (1 << decoder_context->bin_bits[bin])) != 0;
    (*bit) = inv == !b;
    return ICER_RESULT_OK;
}