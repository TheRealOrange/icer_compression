//
// Created by linyi on 19/3/2023.
//

#include "icer.h"

void icer_push_bin_bits(decoder_context_typedef *decoder_context, uint8_t bin, uint16_t bits, uint16_t num_bits) {
    return;
}

int icer_get_bit_from_codeword(decoder_context_typedef *decoder_context, uint8_t bits) {
    uint8_t bitoffset = decoder_context->encode_bit_offset;
    size_t ind = decoder_context->encode_ind;
    uint16_t d, r;
    bitoffset += (bits - 1);
    r = bitoffset / 8;
    d = bitoffset % 8;
    ind += r;

    return (decoder_context->encoded_words[ind] & (1 << bitoffset)) >> bitoffset;
}

int icer_get_bits_from_codeword(decoder_context_typedef *decoder_context, uint8_t bits) {
    int num = 0;
    int bits_to_decode, decoded = 0;
    uint8_t bitoffset = decoder_context->encode_bit_offset;
    size_t ind = decoder_context->encode_ind;
    uint16_t d, r;
    while (bits) {
        bits_to_decode = icer_min_int(bitoffset, bits);
        if (decoder_context->decoded_bits_total + bits_to_decode > decoder_context->encoded_bits_total) {
            return ICER_DECODER_OUT_OF_DATA;
        }
        num |= ((decoder_context->encoded_words[ind] & (((1 << bits_to_decode) - 1) << bitoffset)) >> bitoffset) << decoded;
        decoded += bits_to_decode;
        bitoffset += bits_to_decode;
        r = bitoffset / 8;
        d = bitoffset % 8;
        bitoffset = d;
        if (r) {
            ind++;
        }
    }
    return num;
}

int icer_pop_bits_from_codeword(decoder_context_typedef *decoder_context, uint8_t bits) {
    int num = 0;
    int bits_to_decode, decoded = 0;
    uint16_t d, r;
    while (bits) {
        bits_to_decode = icer_min_int(decoder_context->encode_bit_offset, bits);
        if (decoder_context->decoded_bits_total + bits_to_decode > decoder_context->encoded_bits_total) {
            return ICER_DECODER_OUT_OF_DATA;
        }
        num |= ((decoder_context->encoded_words[decoder_context->encode_ind] & (((1 << bits_to_decode) - 1) << decoder_context->encode_bit_offset)) >> decoder_context->encode_bit_offset) << decoded;
        decoded += bits_to_decode;
        decoder_context->encode_bit_offset += bits_to_decode;
        r = decoder_context->encode_bit_offset / 8;
        d = decoder_context->encode_bit_offset % 8;
        decoder_context->encode_bit_offset = d;
        if (r) {
            decoder_context->encode_ind++;
        }
    }
    return num;
}

int icer_decode_bit(decoder_context_typedef *decoder_context, uint8_t *bit, uint32_t zero_cnt, uint32_t total_cnt) {
    bool inv = false, b;
    uint8_t code_bit;
    uint8_t num_bits = 0;
    uint16_t golomb_k;
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
            icer_get_bit_from_codeword(decoder_context, 1);
            if (code_bit) {

            }
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