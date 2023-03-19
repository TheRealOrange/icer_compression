//
// Created by linyi on 19/3/2023.
//

#include "icer.h"

static inline uint8_t get_bit_category_uint8(const uint8_t* data, uint8_t lsb);
static inline uint8_t get_bit_category_uint16(const uint16_t* data, uint8_t lsb);
static inline bool get_bit_significance_uint8(const uint8_t* data, uint8_t lsb);
static inline bool get_bit_significance_uint16(const uint16_t* data, uint8_t lsb);
static inline int8_t get_sign_uint8(const uint8_t* data, uint8_t lsb);
static inline int8_t get_sign_uint16(const uint16_t* data, uint8_t lsb);

int compress_bitplane_uint8(uint8_t *data, size_t plane_w, size_t plane_h, size_t rowstride,
                            icer_context_model_typedef *context_model,
                            encoder_context_typedef *encoder_context,
                            packet_context* pkt_context) {
    uint8_t *pos;
    uint8_t *rowstart = data;
    int category;
    bool bit;
    uint8_t lsb = pkt_context->lsb;
    uint8_t mask = 0b1 << (lsb-1);

    printf("max out: %zu\n", encoder_context->max_output_length);

    uint8_t *h0, *h1, *v0, *v1, *d0, *d1, *d2, *d3;
    uint8_t h, v, d, tmp;
    int8_t sh0, sh1, sv0, sv1;
    uint8_t sh, sv;
    uint8_t pred_sign, actual_sign, agreement_bit;
    enum icer_pixel_contexts sign_context;
    size_t vert_bound = plane_h - 1;
    size_t hor_bound = plane_w - 1;

    uint8_t prev_plane = lsb+1;
    if (prev_plane >= 8) return ICER_BITPLANE_OUT_OF_RANGE;

    for (size_t row = 0; row < plane_h; row++) {
        pos = rowstart;
        /* keep track of 8 adjacent pixels */
        h0 = pos - 1; h1 = pos + 1;
        v0 = pos - rowstride; v1 = pos + rowstride;
        d0 = v0 - 1; d1 = v0 + 1;
        d2 = v1 - 1; d3 = v1 + 1;

        enum icer_pixel_contexts cntxt;

        for (size_t col = 0; col < plane_w; col++) {
            category = get_bit_category_uint8(pos, lsb);
            bit = ((*pos) & mask) != 0;

            if (category == ICER_CATEGORY_3) {
                /* pass to uncoded bin */
                icer_encode_bit(encoder_context, bit, 1, 2);
            } else {
                if (category == ICER_CATEGORY_0 || category == ICER_CATEGORY_1) {
                    h = 0;
                    v = 0;
                    d = 0;

                    /* consider the horizontally adjacent pixels */
                    if (col > 0) h += get_bit_significance_uint8(h0, lsb) & 1;
                    if (col < hor_bound) h += get_bit_significance_uint8(h1, prev_plane) & 1;

                    /* consider the vertically adjacent pixels */
                    if (row > 0) v += get_bit_significance_uint8(v0, lsb) & 1;
                    if (row < vert_bound) v += get_bit_significance_uint8(v1, prev_plane) & 1;

                    /* consider the diagonally adjacent pixels */
                    if (col > 0 && row > 0) d += get_bit_significance_uint8(d0, lsb) & 1;
                    if (col > 0 && row < vert_bound) d += get_bit_significance_uint8(d2, prev_plane) & 1;
                    if (col < hor_bound && row > 0) d += get_bit_significance_uint8(d1, lsb) & 1;
                    if (col < hor_bound && row < vert_bound) d += get_bit_significance_uint8(d3, prev_plane) & 1;
                }

                if (category == ICER_CATEGORY_0) {
                    if (context_model->subband_type == ICER_SUBBAND_HL) {
                        tmp = h;
                        h = v;
                        v = tmp;
                    }

                    if (context_model->subband_type != ICER_SUBBAND_HH) {
                        cntxt = icer_context_table_ll_lh_hl[h][v][d];
                    } else {
                        cntxt = icer_context_table_hh[h + v][d];
                    }
                } else if (category == ICER_CATEGORY_1) {
                    cntxt = (h + v == 0) ? ICER_CONTEXT_9 : ICER_CONTEXT_10;
                } else if (category == ICER_CATEGORY_2) {
                    cntxt = ICER_CONTEXT_11;
                }

                if (icer_encode_bit(encoder_context, bit, context_model->zero_count[cntxt], context_model->total_count[cntxt]) == ICER_BYTE_QUOTA_EXCEEDED) {
                    return ICER_BYTE_QUOTA_EXCEEDED;
                }

                context_model->total_count[cntxt]++;
                context_model->zero_count[cntxt] += (uint32_t)(!bit);
                if (context_model->total_count[cntxt] >= ICER_CONTEXT_RESCALING_CAP) {
                    context_model->total_count[cntxt] >>= 1;
                    if (context_model->zero_count[cntxt] > context_model->total_count[cntxt]) context_model->zero_count[cntxt] >>= 1;
                    else icer_ceil_div_uint32(context_model->zero_count[cntxt], 2);
                }

                if (category == ICER_CATEGORY_0 && bit) {
                    /* bit is the first magnitude bit to be encoded, thus we will encode the sign bit */
                    /* predict sign bit */
                    sh0 = 0; sh1 = 0;
                    sv0 = 0; sv1 = 0;

                    /* consider the horizontally adjacent pixels */
                    if (col > 0) sh0 = get_sign_uint8(h0, lsb);
                    if (col < hor_bound) sh1 = get_sign_uint8(h1, prev_plane);

                    /* consider the vertically adjacent pixels */
                    if (row > 0) sv0 = get_sign_uint8(v0, lsb);
                    if (row < vert_bound) sv1 = get_sign_uint8(v1, prev_plane);

                    sh = sh0 + sh1 + 2;
                    sv = sv0 + sv1 + 2;

                    if (context_model->subband_type == ICER_SUBBAND_HL) {
                        tmp = sh;
                        sh = sv;
                        sv = tmp;
                    }

                    sign_context = icer_sign_context_table[sh][sv];
                    pred_sign = icer_sign_prediction_table[sh][sv];
                    actual_sign = ((*pos) & 0x80) != 0;

                    agreement_bit = (pred_sign ^ actual_sign) & 1;

                    if (icer_encode_bit(encoder_context, agreement_bit, context_model->zero_count[sign_context], context_model->total_count[sign_context]) == ICER_BYTE_QUOTA_EXCEEDED) {
                        return ICER_BYTE_QUOTA_EXCEEDED;
                    }

                    context_model->total_count[sign_context]++;
                    context_model->zero_count[sign_context] += (uint32_t)(agreement_bit == 0);
                    if (context_model->total_count[sign_context] >= ICER_CONTEXT_RESCALING_CAP) {
                        context_model->total_count[sign_context] >>= 1;
                        if (context_model->zero_count[sign_context] > context_model->total_count[sign_context]) context_model->zero_count[sign_context] >>= 1;
                        else icer_ceil_div_uint32(context_model->zero_count[sign_context], 2);
                    }
                }
            }

            pos++;
            /* increment position pointers of 8 adjacent pixels */
            h0++; h1++; v0++; v1++; d0++; d1++; d2++; d2++;
        }
        rowstart += rowstride;
    }
    if (flush_encode(encoder_context) == ICER_BYTE_QUOTA_EXCEEDED) {
        printf("end flush exceeded output size: %zu bytes\n", encoder_context->output_ind);
        return ICER_BYTE_QUOTA_EXCEEDED;
    }
    return ICER_RESULT_OK;
}

int decompress_bitplane_uint8(uint8_t *data, size_t plane_w, size_t plane_h, size_t rowstride,
                              icer_context_model_typedef *context_model,
                              decoder_context_typedef *decoder_context,
                              packet_context* pkt_context) {
    uint8_t *pos;
    uint8_t *rowstart = data;
    int category;
    uint8_t bit;
    uint8_t lsb = pkt_context->lsb;

    uint8_t *h0, *h1, *v0, *v1, *d0, *d1, *d2, *d3;
    uint8_t h, v, d, tmp;
    int8_t sh0, sh1, sv0, sv1;
    uint8_t sh, sv;
    uint8_t pred_sign, actual_sign, agreement_bit;
    enum icer_pixel_contexts sign_context;
    size_t vert_bound = plane_h - 1;
    size_t hor_bound = plane_w - 1;

    uint8_t prev_plane = lsb+1;
    if (prev_plane >= 8) return ICER_BITPLANE_OUT_OF_RANGE;

    for (size_t row = 0; row < plane_h; row++) {
        pos = rowstart;
        /* keep track of 8 adjacent pixels */
        h0 = pos - 1; h1 = pos + 1;
        v0 = pos - rowstride; v1 = pos + rowstride;
        d0 = v0 - 1; d1 = v0 + 1;
        d2 = v1 - 1; d3 = v1 + 1;

        enum icer_pixel_contexts cntxt;

        for (size_t col = 0; col < plane_w; col++) {
            category = get_bit_category_uint8(pos, lsb);

            if (category == ICER_CATEGORY_3) {
                /* pass to uncoded bin */
                icer_decode_bit(decoder_context, &bit, 1, 2);
                (*pos) |= bit << lsb;
            } else {
                if (category == ICER_CATEGORY_0 || category == ICER_CATEGORY_1) {
                    h = 0;
                    v = 0;
                    d = 0;

                    /* consider the horizontally adjacent pixels */
                    if (col > 0) h += get_bit_significance_uint8(h0, lsb) & 1;
                    if (col < hor_bound) h += get_bit_significance_uint8(h1, prev_plane) & 1;

                    /* consider the vertically adjacent pixels */
                    if (row > 0) v += get_bit_significance_uint8(v0, lsb) & 1;
                    if (row < vert_bound) v += get_bit_significance_uint8(v1, prev_plane) & 1;

                    /* consider the diagonally adjacent pixels */
                    if (col > 0 && row > 0) d += get_bit_significance_uint8(d0, lsb) & 1;
                    if (col > 0 && row < vert_bound) d += get_bit_significance_uint8(d2, prev_plane) & 1;
                    if (col < hor_bound && row > 0) d += get_bit_significance_uint8(d1, lsb) & 1;
                    if (col < hor_bound && row < vert_bound) d += get_bit_significance_uint8(d3, prev_plane) & 1;
                }

                if (category == ICER_CATEGORY_0) {
                    if (context_model->subband_type == ICER_SUBBAND_HL) {
                        tmp = h;
                        h = v;
                        v = tmp;
                    }

                    if (context_model->subband_type != ICER_SUBBAND_HH) {
                        cntxt = icer_context_table_ll_lh_hl[h][v][d];
                    } else {
                        cntxt = icer_context_table_hh[h + v][d];
                    }
                } else if (category == ICER_CATEGORY_1) {
                    cntxt = (h + v == 0) ? ICER_CONTEXT_9 : ICER_CONTEXT_10;
                } else if (category == ICER_CATEGORY_2) {
                    cntxt = ICER_CONTEXT_11;
                }

                icer_decode_bit(decoder_context, &bit, context_model->zero_count[cntxt], context_model->total_count[cntxt]);
                (*pos) |= bit << lsb;

                context_model->total_count[cntxt]++;
                context_model->zero_count[cntxt] += (uint32_t)(!bit);
                if (context_model->total_count[cntxt] >= ICER_CONTEXT_RESCALING_CAP) {
                    context_model->total_count[cntxt] >>= 1;
                    if (context_model->zero_count[cntxt] > context_model->total_count[cntxt]) context_model->zero_count[cntxt] >>= 1;
                    else icer_ceil_div_uint32(context_model->zero_count[cntxt], 2);
                }

                if (category == ICER_CATEGORY_0 && bit) {
                    /* bit is the first magnitude bit to be encoded, thus we will encode the sign bit */
                    /* predict sign bit */
                    sh0 = 0; sh1 = 0;
                    sv0 = 0; sv1 = 0;

                    /* consider the horizontally adjacent pixels */
                    if (col > 0) sh0 = get_sign_uint8(h0, lsb);
                    if (col < hor_bound) sh1 = get_sign_uint8(h1, prev_plane);

                    /* consider the vertically adjacent pixels */
                    if (row > 0) sv0 = get_sign_uint8(v0, lsb);
                    if (row < vert_bound) sv1 = get_sign_uint8(v1, prev_plane);

                    sh = sh0 + sh1 + 2;
                    sv = sv0 + sv1 + 2;

                    if (context_model->subband_type == ICER_SUBBAND_HL) {
                        tmp = sh;
                        sh = sv;
                        sv = tmp;
                    }

                    sign_context = icer_sign_context_table[sh][sv];
                    pred_sign = icer_sign_prediction_table[sh][sv];

                    icer_decode_bit(decoder_context, &agreement_bit, 1, 2);

                    actual_sign = (agreement_bit ^ pred_sign) & 1;
                    (*pos) |= actual_sign << 7;

                    context_model->total_count[sign_context]++;
                    context_model->zero_count[sign_context] += (uint32_t)(agreement_bit == 0);
                    if (context_model->total_count[sign_context] >= ICER_CONTEXT_RESCALING_CAP) {
                        context_model->total_count[sign_context] >>= 1;
                        if (context_model->zero_count[sign_context] > context_model->total_count[sign_context]) context_model->zero_count[sign_context] >>= 1;
                        else icer_ceil_div_uint32(context_model->zero_count[sign_context], 2);
                    }
                }
            }

            pos++;
            /* increment position pointers of 8 adjacent pixels */
            h0++; h1++; v0++; v1++; d0++; d1++; d2++; d2++;
        }
        rowstart += rowstride;
    }
    return ICER_RESULT_OK;
}


void init_context_model_vals(icer_context_model_typedef* context_model, enum icer_subband_types subband_type) {
    context_model->subband_type = subband_type;
    for (size_t i = 0;i <= ICER_CONTEXT_MAX;i++) {
        context_model->zero_count[i] = ICER_DEFAULT_CONTEXT_ZERO_COUNT;
        context_model->total_count[i] = ICER_DEFAULT_CONTEXT_TOTAL_COUNT;
    }
}

static inline uint8_t get_bit_category_uint8(const uint8_t* data, uint8_t lsb) {
    return icer_min_int(__builtin_popcount(((*data) & 0x7f) >> lsb), 3);
}

static inline uint8_t get_bit_category_uint16(const uint16_t* data, uint8_t lsb) {
    return icer_min_int(__builtin_popcount(((*data) & 0x7fff) >> lsb), 3);
}

static inline bool get_bit_significance_uint8(const uint8_t* data, uint8_t lsb) {
    return __builtin_popcount(((*data) & 0x7f) >> lsb) != 0;
}

static inline bool get_bit_significance_uint16(const uint16_t* data, uint8_t lsb) {
    return __builtin_popcount(((*data) & 0x7fff) >> lsb) != 0;
}

static inline int8_t get_sign_uint8(const uint8_t* data, uint8_t lsb) {
    return ((int8_t)(*data) >> 7) * (int8_t)get_bit_significance_uint8(data, lsb);
}

static inline int8_t get_sign_uint16(const uint16_t* data, uint8_t lsb) {
    return ((int16_t)(*data) >> 15) * (int8_t)get_bit_significance_uint16(data, lsb);
}