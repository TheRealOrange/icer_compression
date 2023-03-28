//
// Created by linyi on 19/3/2023.
//

#include "icer.h"

#ifdef USE_ENCODE_FUNCTIONS
#ifndef USER_PROVIDED_BUFFERS
uint16_t icer_encode_circ_buf[ICER_CIRC_BUF_SIZE];
#endif

static inline uint16_t pop_buf(icer_encoder_context_typedef *cntxt);
static inline int16_t alloc_buf(icer_encoder_context_typedef *cntxt);

void icer_init_entropy_coder_context(icer_encoder_context_typedef *encoder_context, uint16_t *encode_buffer, size_t buffer_length, uint8_t *encoder_out, size_t enc_out_max) {
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
    encoder_context->output_buffer[encoder_context->output_ind] = 0;
}


int icer_encode_bit(icer_encoder_context_typedef *encoder_context, uint8_t bit, uint32_t zero_cnt, uint32_t total_cnt) {

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
    uint16_t golomb_out_code;
    uint8_t golomb_out_bits;
    uint16_t bit16 = (bit != 0);

    if (encoder_context->bin_current_buf[bin] == -1) {
        encoder_context->bin_current_buf[bin] = alloc_buf(encoder_context);
        if (encoder_context->bin_current_buf[bin] == -1) {
            if (icer_flush_encode(encoder_context) == ICER_BYTE_QUOTA_EXCEEDED) return ICER_BYTE_QUOTA_EXCEEDED;
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
            golomb_out_code = golomb_k + ((golomb_k < icer_golomb_coders[bin].i) ? 0 : icer_golomb_coders[bin].i);
            golomb_out_bits = icer_golomb_coders[bin].l + (golomb_k >= icer_golomb_coders[bin].i);
            icer_reverse_bits(&golomb_out_code, golomb_out_bits);
            (*curr_bin) = (golomb_out_bits << ICER_ENC_BUF_BITS_OFFSET);
            (*curr_bin) |= (golomb_out_code & ICER_ENC_BUF_DATA_MASK);
            (*curr_bin) |= ICER_ENC_BUF_DONE_MASK;
            encoder_context->bin_current_buf[bin] = -1;
        } else if (((*curr_bin) & ICER_ENC_BUF_DATA_MASK) >= icer_golomb_coders[bin].m) {
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
        if (icer_custom_coding_scheme[bin][prefix].input_code_bits == encoder_context->bin_current_buf_bits[bin]) {
            (*curr_bin) = icer_custom_coding_scheme[bin][prefix].output_code_bits << ICER_ENC_BUF_BITS_OFFSET;
            (*curr_bin) |= icer_custom_coding_scheme[bin][prefix].output_code & ICER_ENC_BUF_DATA_MASK;
            (*curr_bin) |= ICER_ENC_BUF_DONE_MASK;
            encoder_context->bin_current_buf[bin] = -1;
            encoder_context->bin_current_buf_bits[bin] = 0;
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

int icer_popbuf_while_avail(icer_encoder_context_typedef *encoder_context) {
    uint16_t out, bits;
    uint16_t d, r;
    int bits_to_encode;
    while (encoder_context->used > 0 && (encoder_context->encode_buffer[encoder_context->head] & ICER_ENC_BUF_DONE_MASK)) {
        out = pop_buf(encoder_context);
        bits = out >> ICER_ENC_BUF_BITS_OFFSET;
        while (bits) {
            bits_to_encode = icer_min_int(8-encoder_context->output_bit_offset, bits);
            encoder_context->output_buffer[encoder_context->output_ind] |= (out & ((1 << bits_to_encode) - 1)) << (encoder_context->output_bit_offset);
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

int icer_flush_encode(icer_encoder_context_typedef *encoder_context) {
    uint16_t *first = encoder_context->encode_buffer + encoder_context->head;
    icer_custom_flush_typedef *flush;
    uint16_t prefix;
    if (((*first) & ICER_ENC_BUF_DONE_MASK) == 0) {
        uint8_t bin = (*first) >> ICER_ENC_BUF_BITS_OFFSET;

        uint16_t golomb_k;
        uint16_t golomb_out_code;
        uint8_t golomb_out_bits;

        if (bin > ICER_ENC_BIN_8) {
            /* golomb code bins */
            golomb_k = (*first) & ICER_ENC_BUF_DATA_MASK;
            if (golomb_k == icer_golomb_coders[bin].m - 1) {
                *first = 1 << ICER_ENC_BUF_BITS_OFFSET;
                *first |= 1;
                *first |= ICER_ENC_BUF_DONE_MASK;
            } else {
                golomb_out_code = golomb_k + ((golomb_k < icer_golomb_coders[bin].i) ? 0 : icer_golomb_coders[bin].i);
                golomb_out_bits = icer_golomb_coders[bin].l + (golomb_k >= icer_golomb_coders[bin].i);
                icer_reverse_bits(&golomb_out_code, golomb_out_bits);
                *first = (golomb_out_bits << ICER_ENC_BUF_BITS_OFFSET);
                *first |= (golomb_out_code & ICER_ENC_BUF_DATA_MASK);
                *first |= ICER_ENC_BUF_DONE_MASK;
            }
            encoder_context->bin_current_buf[bin] = -1;
        } else if (bin != ICER_ENC_BIN_1) {
            /* custom non prefix code bins */
            flush = &(icer_custom_code_flush_bits[bin][(*first) & ICER_ENC_BUF_DATA_MASK][encoder_context->bin_current_buf_bits[bin]]);
            (*first) |= flush->flush_bit << encoder_context->bin_current_buf_bits[bin];
            encoder_context->bin_current_buf_bits[bin] += flush->flush_bit_numbers;

            prefix = (*first) & ICER_ENC_BUF_DATA_MASK;

            (*first) = icer_custom_coding_scheme[bin][prefix].output_code_bits << ICER_ENC_BUF_BITS_OFFSET;
            (*first) |= icer_custom_coding_scheme[bin][prefix].output_code & ICER_ENC_BUF_DATA_MASK;
            (*first) |= ICER_ENC_BUF_DONE_MASK;
            encoder_context->bin_current_buf[bin] = -1;
            encoder_context->bin_current_buf_bits[bin] = 0;
        } else {
            /* uncoded bin */
            // this should never happen?
        }
    }

    if (icer_popbuf_while_avail(encoder_context) == ICER_BYTE_QUOTA_EXCEEDED) return ICER_BYTE_QUOTA_EXCEEDED;
    return ICER_RESULT_OK;
}

/* circular buffer helper functions */

static inline uint16_t pop_buf(icer_encoder_context_typedef *cntxt) {
    if (cntxt->used > 0) cntxt->used--;
    uint16_t res = cntxt->encode_buffer[cntxt->head];
    cntxt->head = (cntxt->head + 1) % cntxt->buffer_length;
    return res;
}

static inline int16_t alloc_buf(icer_encoder_context_typedef *cntxt) {
    if (cntxt->used >= cntxt->buffer_length) return -1;
    cntxt->used++;
    int16_t ind = (int16_t)cntxt->tail;
    cntxt->tail = (cntxt->tail+1) % cntxt->buffer_length;
    return ind;
}

/* data packet functions */

int icer_allocate_data_packet(icer_image_segment_typedef **pkt, icer_output_data_buf_typedef * const output_data, uint8_t segment_num, const icer_packet_context *context) {
    size_t buf_len = output_data->size_allocated - output_data->size_used;
    if (buf_len < sizeof(icer_image_segment_typedef)) {
        return ICER_BYTE_QUOTA_EXCEEDED;
    }
    (*pkt) = (icer_image_segment_typedef *) (output_data->data_start + output_data->size_used);
    (*pkt)->preamble = ICER_PACKET_PREAMBLE;
    (*pkt)->decomp_level = context->decomp_level;
    (*pkt)->subband_type = context->subband_type;
    (*pkt)->segment_number = segment_num;
    (*pkt)->lsb_chan = context->lsb | ICER_SET_CHANNEL_MACRO(context->channel);
    (*pkt)->ll_mean_val = context->ll_mean_val;
    (*pkt)->image_w = context->image_w;
    (*pkt)->image_h = context->image_h;
    (*pkt)->data_crc32 = 0;
    (*pkt)->crc32 = 0;

    output_data->size_used += sizeof(icer_image_segment_typedef);
    buf_len -= sizeof(icer_image_segment_typedef);

    // store max data length first
    (*pkt)->data_length = buf_len;

    return ICER_RESULT_OK;
}

#endif