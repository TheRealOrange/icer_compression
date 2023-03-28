//
// Created by linyi on 19/3/2023.
//

#include <string.h>
#include "icer.h"

#ifdef USE_DECODE_FUNCTIONS

#define ICER_BITMASK_MACRO(x) (((unsigned)1 << x) - 1)

void icer_init_entropy_decoder_context(icer_decoder_context_typedef *decoder_context, uint8_t *encoded_words, size_t encoded_bits) {
    decoder_context->encoded_bits_total = encoded_bits;
    decoder_context->decoded_bits_total = 0;
    decoder_context->decoded_words = 0;

    decoder_context->encoded_words = encoded_words;
    decoder_context->encode_bit_offset = 0;
    decoder_context->encode_ind = 0;

    for (size_t it = 0;it < ICER_ENCODER_BIN_MAX+1;it++) {
        decoder_context->bin_bits[it] = 0;
        decoder_context->bin_decode_index[it] = 0;
        for (size_t j = 0;j < ICER_DECODER_BIT_BIN_MAX;j++) decoder_context->bin_buf[it][j] = 0;
    }
}

void icer_push_bin_bits(icer_decoder_context_typedef *decoder_context, uint8_t bin, uint16_t bits, uint16_t num_bits) {
    int32_t bin_ind, bin_bit_offset;
    bin_ind = decoder_context->bin_bits[bin] / 32;
    bin_bit_offset = decoder_context->bin_bits[bin] % 32;
    decoder_context->bin_bits[bin] += num_bits;

    int bits_to_push;
    while (num_bits) {
        bits_to_push = icer_min_int(num_bits, 32-bin_bit_offset);
        decoder_context->bin_buf[bin][bin_ind] |= ((bits & ICER_BITMASK_MACRO(bits_to_push)) << bin_bit_offset);
        num_bits -= bits_to_push;

        bin_bit_offset += bits_to_push;
        bin_ind += bin_bit_offset / 32;
        bin_bit_offset %= 32;
    }
}

int icer_get_bit_from_codeword(icer_decoder_context_typedef *decoder_context, uint8_t bits) {
    uint8_t bitoffset = decoder_context->encode_bit_offset;
    size_t ind = decoder_context->encode_ind;
    uint16_t d, r;
    bitoffset += (bits - 1);
    r = bitoffset / 8;
    d = bitoffset % 8;
    ind += r;
    bitoffset = d;

    return (decoder_context->encoded_words[ind] & (1 << bitoffset)) >> bitoffset;
}

int icer_get_bits_from_codeword(icer_decoder_context_typedef *decoder_context, uint8_t bits) {
    int num = 0;
    int bits_to_decode, decoded = 0;
    uint8_t bitoffset = decoder_context->encode_bit_offset;
    size_t ind = decoder_context->encode_ind;
    uint16_t d, r;
    while (bits) {
        bits_to_decode = icer_min_int(8-bitoffset, bits);
        if (decoder_context->decoded_bits_total + bits_to_decode > decoder_context->encoded_bits_total) {
            return ICER_DECODER_OUT_OF_DATA;
        }
        num |= (int)((((decoder_context->encoded_words[ind] & (ICER_BITMASK_MACRO(bits_to_decode)) << bitoffset)) >> bitoffset) << decoded);
        bits -= bits_to_decode;
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

int icer_pop_bits_from_codeword(icer_decoder_context_typedef *decoder_context, uint8_t bits) {
    int num = 0;
    int bits_to_decode, decoded = 0;
    uint16_t d, r;
    while (bits) {
        bits_to_decode = icer_min_int(8-decoder_context->encode_bit_offset, bits);
        if (decoder_context->decoded_bits_total + bits_to_decode > decoder_context->encoded_bits_total) {
            return ICER_DECODER_OUT_OF_DATA;
        }
        num |= (int)(((decoder_context->encoded_words[decoder_context->encode_ind] & (ICER_BITMASK_MACRO(bits_to_decode) << decoder_context->encode_bit_offset)) >> decoder_context->encode_bit_offset) << decoded);
        bits -= bits_to_decode;
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


int icer_decode_bit(icer_decoder_context_typedef *decoder_context, uint8_t *bit, uint32_t zero_cnt, uint32_t total_cnt) {
    bool inv = false, b;
    int code_bit;
    uint16_t codeword;
    uint8_t num_bits;
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

    if (decoder_context->bin_bits[bin] <= 0 || decoder_context->decoded_words - decoder_context->bin_decode_index[bin] >= ICER_CIRC_BUF_SIZE) {
        /* ran out of bits in the bit, time to process a new codeword */
        decoder_context->bin_bits[bin] = 0;
        memset(decoder_context->bin_buf[bin], 0, ICER_DECODER_BIT_BIN_MAX*sizeof(decoder_context->bin_buf[bin][0]));
        if (bin > ICER_ENC_BIN_8) {
            /* golomb code bins */
            code_bit = icer_get_bit_from_codeword(decoder_context, 1);
            if (code_bit == ICER_DECODER_OUT_OF_DATA) return ICER_DECODER_OUT_OF_DATA; //major oops moment
            if (code_bit) {
                /* if the first bit is one, return m 0s */
                icer_pop_bits_from_codeword(decoder_context, 1);
                icer_push_bin_bits(decoder_context, bin, 0b0, icer_golomb_coders[bin].m);
            } else {
                golomb_k = icer_get_bits_from_codeword(decoder_context, icer_golomb_coders[bin].l);
                icer_reverse_bits(&golomb_k, icer_golomb_coders[bin].l);
                if (golomb_k < icer_golomb_coders[bin].i) {
                    icer_pop_bits_from_codeword(decoder_context, icer_golomb_coders[bin].l);
                    icer_push_bin_bits(decoder_context, bin, 0b1, 1);
                    icer_push_bin_bits(decoder_context, bin, 0b0, golomb_k);
                } else {
                    golomb_k = icer_pop_bits_from_codeword(decoder_context, icer_golomb_coders[bin].l + 1);
                    icer_reverse_bits(&golomb_k, icer_golomb_coders[bin].l + 1);
                    icer_push_bin_bits(decoder_context, bin, 0b1, 1);
                    icer_push_bin_bits(decoder_context, bin, 0b0, golomb_k - icer_golomb_coders[bin].i);
                }
            }
        } else if (bin != ICER_ENC_BIN_1) {
            /* custom non prefix code bins */
            codeword = 0;
            num_bits = 0;
            do {
                if (decoder_context->decoded_bits_total + num_bits + 1 >= decoder_context->encoded_bits_total) return ICER_DECODER_OUT_OF_DATA;
                codeword |= icer_get_bit_from_codeword(decoder_context, num_bits+1) << num_bits;
                num_bits++;
                if (codeword < 32) {
                    if (icer_custom_decode_scheme[bin][codeword].input_code_bits == num_bits) {
                        icer_push_bin_bits(decoder_context, bin, icer_custom_decode_scheme[bin][codeword].output_code, icer_custom_decode_scheme[bin][codeword].output_code_bits);
                        int test = icer_pop_bits_from_codeword(decoder_context, num_bits);
                        if (codeword != test) {
                            return ICER_DECODED_INVALID_DATA;
                        } else {
                            break;
                        }
                    }
                } else {
                    return ICER_DECODED_INVALID_DATA;
                }
            } while (num_bits < 10);
        } else {
            /* uncoded bin */
            code_bit = icer_pop_bits_from_codeword(decoder_context, 1);
            if (code_bit == ICER_DECODER_OUT_OF_DATA) return ICER_DECODER_OUT_OF_DATA;
            icer_push_bin_bits(decoder_context, bin, code_bit != 0, 1);
        }

        decoder_context->decoded_words++;
        decoder_context->bin_decode_index[bin] = decoder_context->decoded_words;
    }
    int32_t bin_ind = decoder_context->bin_bits[bin] / 32;
    int32_t bit_offset = decoder_context->bin_bits[bin] % 32;
    b = (decoder_context->bin_buf[bin][bin_ind] & (1 << (bit_offset-1))) != 0;
    decoder_context->bin_buf[bin][bin_ind] &= ~(1 << (bit_offset-1));
    decoder_context->bin_bits[bin]--;
    (*bit) = inv == !b;

    return ICER_RESULT_OK;
}

#endif