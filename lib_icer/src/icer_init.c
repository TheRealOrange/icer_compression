//
// Created by linyi on 19/3/2023.
//

#include "icer.h"

#define INIT_CODING_SCHEME(bin, inp, inp_bits, out, out_bits) { \
icer_custom_coding_scheme[bin][inp].input_code_bits = inp_bits;       \
icer_custom_coding_scheme[bin][inp].output_code = out;                \
icer_custom_coding_scheme[bin][inp].output_code_bits = out_bits;      \
}

#define INIT_FLUSH_BITS(bin, inp, inp_bits, out, out_bits) { \
icer_custom_code_flush_bits[bin][inp][inp_bits].flush_bit = out;   \
icer_custom_code_flush_bits[bin][inp][inp_bits].flush_bit_numbers = out_bits; \
}

#define INIT_DECODE_SCHEME(bin, out, out_bits, inp, inp_bits) { \
icer_custom_decode_scheme[bin][inp].input_code_bits = inp_bits;       \
icer_custom_decode_scheme[bin][inp].output_code = out;                \
icer_custom_decode_scheme[bin][inp].output_code_bits = out_bits;      \
}

int icer_init() {
    icer_init_golombcoder();
#ifdef USE_ENCODE_FUNCTIONS
    icer_init_codingscheme();
#endif
#ifdef USE_DECODE_FUNCTIONS
    icer_init_decodescheme();
#endif
    icer_init_flushbits();

    return ICER_RESULT_OK;
}

#ifdef USE_DECODE_FUNCTIONS
void icer_init_decodescheme() {
    for (int it = 0; it <= ICER_ENCODER_BIN_MAX; it++) {
        for (int j = 0; j < CUSTOM_CODING_MAX_LOOKUP; j++) {
            icer_custom_decode_scheme[it][j].input_code_bits = 0;
            icer_custom_decode_scheme[it][j].output_code_bits = 0;
        }
    }

    INIT_DECODE_SCHEME(ICER_ENC_BIN_2, 0b01, 2, 0b10, 2);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_2, 0b011, 3, 0b011, 3);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_2, 0b0111, 4, 0b1111, 4);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_2, 0b1111, 4, 0b10000, 5);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_2, 0b10, 2, 0b01, 2);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_2, 0b100, 3, 0b100, 3);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_2, 0b1000, 4, 0b1000, 4);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_2, 0b10000, 5, 0b00000, 5);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_2, 0b00000, 5, 0b0111, 4);

    INIT_DECODE_SCHEME(ICER_ENC_BIN_3, 0b10, 2, 0b01, 2);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_3, 0b100, 3, 0b00, 2);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_3, 0b0000, 4, 0b011, 3);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_3, 0b11000, 5, 0b10010, 5);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_3, 0b01000, 5, 0b1111, 4);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_3, 0b01, 2, 0b110, 3);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_3, 0b0011, 4, 0b0111, 4);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_3, 0b1011, 4, 0b00010, 5);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_3, 0b111, 3, 0b1010, 4);

    INIT_DECODE_SCHEME(ICER_ENC_BIN_4, 0b10, 2, 0b10, 2);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_4, 0b100, 3, 0b011, 3);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_4, 0b000, 3, 0b00, 2);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_4, 0b01, 2, 0b01, 2);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_4, 0b11, 2, 0b111, 3);

    INIT_DECODE_SCHEME(ICER_ENC_BIN_5, 0b00, 2, 0b1, 1);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_5, 0b010, 3, 0b000, 3);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_5, 0b110, 3, 0b1010, 4);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_5, 0b101, 3, 0b0010, 4);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_5, 0b1001, 4, 0b1110, 4);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_5, 0b00001, 5, 0b0100, 4);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_5, 0b10001, 5, 0b00110, 5);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_5, 0b011, 3, 0b1100, 4);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_5, 0b111, 3, 0b10110, 5);

    INIT_DECODE_SCHEME(ICER_ENC_BIN_6, 0b1, 1, 0b10, 2);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_6, 0b010, 3, 0b011, 3);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_6, 0b110, 3, 0b1111, 4);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_6, 0b100, 3, 0b101, 3);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_6, 0b1000, 4, 0b001, 3);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_6, 0b10000, 5, 0b0111, 4);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_6, 0b00000, 5, 0b00, 2);

    INIT_DECODE_SCHEME(ICER_ENC_BIN_7, 0b000, 3, 0b0, 1);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_7, 0b100, 3, 0b001, 3);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_7, 0b010, 3, 0b101, 3);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_7, 0b110, 3, 0b01111, 5);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_7, 0b11, 2, 0b0111, 4);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_7, 0b001, 3, 0b011, 3);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_7, 0b101, 3, 0b11111, 5);

    INIT_DECODE_SCHEME(ICER_ENC_BIN_8, 0b10, 2, 0b101, 3);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_8, 0b100, 3, 0b001, 3);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_8, 0b0000, 4, 0b0, 1);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_8, 0b01000, 5, 0b0111, 4);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_8, 0b11000, 5, 0b01111, 5);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_8, 0b01, 2, 0b011, 3);
    INIT_DECODE_SCHEME(ICER_ENC_BIN_8, 0b11, 2, 0b11111, 5);


    for (int it = 0; it <= ICER_ENCODER_BIN_MAX; it++) {
        for (int j = 0; j < CUSTOM_CODING_MAX_LOOKUP; j++) {
            if (icer_custom_decode_scheme[it][j].output_code_bits != 0) {
                uint8_t reversed = 0;
                for (int b = 0; b < icer_custom_decode_scheme[it][j].output_code_bits; b++) {
                    reversed <<= 1;
                    reversed |= icer_custom_decode_scheme[it][j].output_code & 1;
                    icer_custom_decode_scheme[it][j].output_code >>= 1;
                }
                icer_custom_decode_scheme[it][j].output_code = reversed;
            }
        }
    }
}
#endif

#ifdef USE_ENCODE_FUNCTIONS
void icer_init_codingscheme() {
    for (int it = 0; it <= ICER_ENCODER_BIN_MAX; it++) {
        for (int j = 0; j < CUSTOM_CODING_MAX_LOOKUP; j++) icer_custom_coding_scheme[it][j].input_code_bits = 0;
    }

    INIT_CODING_SCHEME(ICER_ENC_BIN_2, 0b01, 2, 0b10, 2);
    INIT_CODING_SCHEME(ICER_ENC_BIN_2, 0b011, 3, 0b011, 3);
    INIT_CODING_SCHEME(ICER_ENC_BIN_2, 0b0111, 4, 0b1111, 4);
    INIT_CODING_SCHEME(ICER_ENC_BIN_2, 0b1111, 4, 0b10000, 5);
    INIT_CODING_SCHEME(ICER_ENC_BIN_2, 0b10, 2, 0b01, 2);
    INIT_CODING_SCHEME(ICER_ENC_BIN_2, 0b100, 3, 0b100, 3);
    INIT_CODING_SCHEME(ICER_ENC_BIN_2, 0b1000, 4, 0b1000, 4);
    INIT_CODING_SCHEME(ICER_ENC_BIN_2, 0b10000, 5, 0b00000, 5);
    INIT_CODING_SCHEME(ICER_ENC_BIN_2, 0b00000, 5, 0b0111, 4);

    INIT_CODING_SCHEME(ICER_ENC_BIN_3, 0b10, 2, 0b01, 2);
    INIT_CODING_SCHEME(ICER_ENC_BIN_3, 0b100, 3, 0b00, 2);
    INIT_CODING_SCHEME(ICER_ENC_BIN_3, 0b0000, 4, 0b011, 3);
    INIT_CODING_SCHEME(ICER_ENC_BIN_3, 0b11000, 5, 0b10010, 5);
    INIT_CODING_SCHEME(ICER_ENC_BIN_3, 0b01000, 5, 0b1111, 4);
    INIT_CODING_SCHEME(ICER_ENC_BIN_3, 0b01, 2, 0b110, 3);
    INIT_CODING_SCHEME(ICER_ENC_BIN_3, 0b0011, 4, 0b0111, 4);
    INIT_CODING_SCHEME(ICER_ENC_BIN_3, 0b1011, 4, 0b00010, 5);
    INIT_CODING_SCHEME(ICER_ENC_BIN_3, 0b111, 3, 0b1010, 4);

    INIT_CODING_SCHEME(ICER_ENC_BIN_4, 0b10, 2, 0b10, 2);
    INIT_CODING_SCHEME(ICER_ENC_BIN_4, 0b100, 3, 0b011, 3);
    INIT_CODING_SCHEME(ICER_ENC_BIN_4, 0b000, 3, 0b00, 2);
    INIT_CODING_SCHEME(ICER_ENC_BIN_4, 0b01, 2, 0b01, 2);
    INIT_CODING_SCHEME(ICER_ENC_BIN_4, 0b11, 2, 0b111, 3);

    INIT_CODING_SCHEME(ICER_ENC_BIN_5, 0b00, 2, 0b1, 1);
    INIT_CODING_SCHEME(ICER_ENC_BIN_5, 0b010, 3, 0b000, 3);
    INIT_CODING_SCHEME(ICER_ENC_BIN_5, 0b110, 3, 0b1010, 4);
    INIT_CODING_SCHEME(ICER_ENC_BIN_5, 0b101, 3, 0b0010, 4);
    INIT_CODING_SCHEME(ICER_ENC_BIN_5, 0b1001, 4, 0b1110, 4);
    INIT_CODING_SCHEME(ICER_ENC_BIN_5, 0b00001, 5, 0b0100, 4);
    INIT_CODING_SCHEME(ICER_ENC_BIN_5, 0b10001, 5, 0b00110, 5);
    INIT_CODING_SCHEME(ICER_ENC_BIN_5, 0b011, 3, 0b1100, 4);
    INIT_CODING_SCHEME(ICER_ENC_BIN_5, 0b111, 3, 0b10110, 5);

    INIT_CODING_SCHEME(ICER_ENC_BIN_6, 0b1, 1, 0b10, 2);
    INIT_CODING_SCHEME(ICER_ENC_BIN_6, 0b010, 3, 0b011, 3);
    INIT_CODING_SCHEME(ICER_ENC_BIN_6, 0b110, 3, 0b1111, 4);
    INIT_CODING_SCHEME(ICER_ENC_BIN_6, 0b100, 3, 0b101, 3);
    INIT_CODING_SCHEME(ICER_ENC_BIN_6, 0b1000, 4, 0b001, 3);
    INIT_CODING_SCHEME(ICER_ENC_BIN_6, 0b10000, 5, 0b0111, 4);
    INIT_CODING_SCHEME(ICER_ENC_BIN_6, 0b00000, 5, 0b00, 2);

    INIT_CODING_SCHEME(ICER_ENC_BIN_7, 0b000, 3, 0b0, 1);
    INIT_CODING_SCHEME(ICER_ENC_BIN_7, 0b100, 3, 0b001, 3);
    INIT_CODING_SCHEME(ICER_ENC_BIN_7, 0b010, 3, 0b101, 3);
    INIT_CODING_SCHEME(ICER_ENC_BIN_7, 0b110, 3, 0b01111, 5);
    INIT_CODING_SCHEME(ICER_ENC_BIN_7, 0b11, 2, 0b0111, 4);
    INIT_CODING_SCHEME(ICER_ENC_BIN_7, 0b001, 3, 0b011, 3);
    INIT_CODING_SCHEME(ICER_ENC_BIN_7, 0b101, 3, 0b11111, 5);

    INIT_CODING_SCHEME(ICER_ENC_BIN_8, 0b10, 2, 0b101, 3);
    INIT_CODING_SCHEME(ICER_ENC_BIN_8, 0b100, 3, 0b001, 3);
    INIT_CODING_SCHEME(ICER_ENC_BIN_8, 0b0000, 4, 0b0, 1);
    INIT_CODING_SCHEME(ICER_ENC_BIN_8, 0b01000, 5, 0b0111, 4);
    INIT_CODING_SCHEME(ICER_ENC_BIN_8, 0b11000, 5, 0b01111, 5);
    INIT_CODING_SCHEME(ICER_ENC_BIN_8, 0b01, 2, 0b011, 3);
    INIT_CODING_SCHEME(ICER_ENC_BIN_8, 0b11, 2, 0b11111, 5);
}
#endif

void icer_init_flushbits() {
    INIT_FLUSH_BITS(ICER_ENC_BIN_2, 0b1, 1, 0, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_2, 0b11, 2, 0, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_2, 0b111, 3, 0, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_2, 0b0, 1, 1, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_2, 0b00, 2, 1, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_2, 0b000, 3, 1, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_2, 0b0000, 4, 0, 1);

    INIT_FLUSH_BITS(ICER_ENC_BIN_3, 0b0, 1, 1, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_3, 0b00, 2, 1, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_3, 0b000, 3, 0, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_3, 0b1000, 4, 0, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_3, 0b1, 1, 0, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_3, 0b11, 2, 1, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_3, 0b011, 3, 0, 1);

    INIT_FLUSH_BITS(ICER_ENC_BIN_4, 0b0, 1, 1, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_4, 0b00, 2, 0, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_4, 0b1, 1, 0, 1);

    INIT_FLUSH_BITS(ICER_ENC_BIN_5, 0b0, 1, 0, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_5, 0b10, 2, 0, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_5, 0b01, 2, 1, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_5, 0b001, 3, 1, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_5, 0b0001, 4, 0, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_5, 0b1, 1, 0b01, 2);
    INIT_FLUSH_BITS(ICER_ENC_BIN_5, 0b11, 2, 0, 1);

    INIT_FLUSH_BITS(ICER_ENC_BIN_6, 0b0, 1, 0b01, 2);
    INIT_FLUSH_BITS(ICER_ENC_BIN_6, 0b01, 2, 0, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_6, 0b00, 2, 1, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_6, 0b000, 3, 1, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_6, 0b0000, 4, 0, 1);

    INIT_FLUSH_BITS(ICER_ENC_BIN_7, 0b0, 1, 0b00, 2);
    INIT_FLUSH_BITS(ICER_ENC_BIN_7, 0b00, 2, 0, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_7, 0b10, 2, 0, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_7, 0b1, 1, 1, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_7, 0b01, 2, 0, 1);

    INIT_FLUSH_BITS(ICER_ENC_BIN_8, 0b0, 1, 1, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_8, 0b00, 2, 1, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_8, 0b000, 3, 0, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_8, 0b1000, 4, 0, 1);
    INIT_FLUSH_BITS(ICER_ENC_BIN_8, 0b1, 1, 0, 1);
}

void icer_init_golombcoder() {
    for (int it = 0; it <= ICER_ENCODER_BIN_MAX; it++) {
        if (icer_bin_coding_scheme[it] > 0) {
            unsigned int m = icer_bin_coding_scheme[it];
            icer_golomb_coders[it].m = m;

            // compute ceil( log2( m ) )
            unsigned int l = 31 - __builtin_clz(m);
            l += ((m ^ (1 << l)) != 0);

            // compute 2^l - m
            unsigned int i = icer_pow_uint(2, l) - m;

            icer_golomb_coders[it].i = i;
            icer_golomb_coders[it].l = l;
        }
    }
}