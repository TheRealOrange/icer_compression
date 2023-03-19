//
// Created by linyi on 19/3/2023.
//

#include "icer.h"

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