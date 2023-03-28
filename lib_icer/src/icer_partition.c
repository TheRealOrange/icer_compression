//
// Created by linyi on 19/3/2023.
//

#include "icer.h"

int icer_generate_partition_parameters(partition_param_typdef *params, size_t ll_w, size_t ll_h, uint16_t segments) {
    uint16_t r, c, r_t, h_t, x_t, c_t0, y_t, r_t0;
    uint16_t x_b = 0, c_b0 = 0, y_b = 0, r_b0 = 0;

    if (segments > (ll_w * ll_h) || segments > ICER_MAX_SEGMENTS) {
        return ICER_TOO_MANY_SEGMENTS;
    }

    if (ll_h > (segments - 1) * ll_w) {
        r = segments;
    } else {
        for (r = 1; r < segments && (r + 1) * r * ll_w < ll_h * segments; r++);
    }
    c = segments / r;
    r_t = (c + 1) * r - segments;
    h_t = icer_max_uint(r_t, ((2 * ll_h * c * r_t + segments) / 2) / segments);
    x_t = ll_w / c;
    c_t0 = (x_t + 1) * c - ll_w;
    y_t = h_t / r_t;
    r_t0 = (y_t + 1) * r_t - h_t;

    if (r_t < r) {
        x_b = ll_w / (c + 1);
        c_b0 = (x_b + 1) * (c + 1) - ll_w;
        y_b = (ll_h - h_t) / (r - r_t);
        r_b0 = (y_b + 1) * (r - r_t) - (ll_h - h_t);
    }

    params->w = ll_w;
    params->h = ll_h;
    params->s = segments;

    params->r = r;
    params->c = c;
    params->r_t = r_t;
    params->h_t = h_t;
    params->x_t = x_t;
    params->c_t0 = c_t0;
    params->y_t = y_t;
    params->r_t0 = r_t0;

    params->x_b = x_b;
    params->c_b0 = c_b0;
    params->y_b = y_b;
    params->r_b0 = r_b0;

    return ICER_RESULT_OK;
}

#ifdef USE_UINT8_FUNCTIONS

#ifdef USE_ENCODE_FUNCTIONS
int icer_compress_partition_uint8(const uint8_t *data, partition_param_typdef *params, size_t rowstride, icer_packet_context *pkt_context,
                                  icer_output_data_buf_typedef * const output_data,icer_image_segment_typedef *segments_encoded[]) {
    int res;
    size_t segment_w, segment_h;
    const uint8_t *segment_start;
    uint16_t segment_num = 0;

    size_t partition_col_ind;
    size_t partition_row_ind = 0;

    icer_context_model_typedef context_model;
    icer_encoder_context_typedef context;
    icer_image_segment_typedef *seg;

    uint32_t data_in_bytes;
    /*
     * process top region which consists of c columns
     * height of top region is h_t and it contains r_t rows
     */
    for (uint16_t row = 0; row < params->r_t; row++) {
        /*
         * the first r_t0 rows have height y_t
         * the remainder have height y_t + 1
         */
        segment_h = params->y_t + ((row >= params->r_t0) ? 1 : 0);
        partition_col_ind = 0;

        for (uint16_t col = 0; col < params->c; col++) {
            /* the first c_t0 columns have width x_t
             * the remainder have width x_t + 1
             */
            segment_w = params->x_t + ((col >= params->c_t0) ? 1 : 0);
            segment_start = data + partition_row_ind * rowstride + partition_col_ind;
            partition_col_ind += segment_w;

            icer_init_context_model_vals(&context_model, pkt_context->subband_type);
            res = icer_allocate_data_packet(&seg, output_data, segment_num, pkt_context);
            if (res != ICER_RESULT_OK) return res;

            icer_init_entropy_coder_context(&context, icer_encode_circ_buf, ICER_CIRC_BUF_SIZE,
                                            (uint8_t *) seg + sizeof(icer_image_segment_typedef), seg->data_length);
            res = icer_compress_bitplane_uint8(segment_start, segment_w, segment_h, rowstride, &context_model, &context,
                                             pkt_context);
            if (res != ICER_RESULT_OK) {
                output_data->size_used -= sizeof(icer_image_segment_typedef);
                return res;
            }

            data_in_bytes = context.output_ind + (context.output_bit_offset > 0);
            seg->data_length = context.output_ind * 8 + context.output_bit_offset;
            seg->data_crc32 = icer_calculate_segment_crc32(seg);
            seg->crc32 = icer_calculate_packet_crc32(seg);
            output_data->size_used += data_in_bytes;

            segments_encoded[segment_num] = seg;

            segment_num++;
        }
        partition_row_ind += segment_h;
    }

    /*
     * if the bottom region exists, process bottom region
     * which consists of c+1 columns
     */
    for (uint16_t row = 0; row < (params->r - params->r_t); row++) {
        /*
         * the first r_b0 rows have height y_b
         * the remainder have height y_b + 1
         */
        segment_h = params->y_b + ((row >= params->r_b0) ? 1 : 0);
        partition_col_ind = 0;

        for (uint16_t col = 0; col < (params->c + 1); col++) {
            /* the first c_b0 columns have width x_b
             * the remainder have width x_b + 1
             */
            segment_w = params->x_b + ((col >= params->c_b0) ? 1 : 0);
            segment_start = data + partition_row_ind * rowstride + partition_col_ind;
            partition_col_ind += segment_w;

            icer_init_context_model_vals(&context_model, pkt_context->subband_type);
            res = icer_allocate_data_packet(&seg, output_data, segment_num, pkt_context);
            if (res != ICER_RESULT_OK) return res;

            icer_init_entropy_coder_context(&context, icer_encode_circ_buf, ICER_CIRC_BUF_SIZE,
                                            (uint8_t *) seg + sizeof(icer_image_segment_typedef), seg->data_length);
            res = icer_compress_bitplane_uint8(segment_start, segment_w, segment_h, rowstride, &context_model, &context,
                                             pkt_context);
            if (res != ICER_RESULT_OK) {
                output_data->size_used -= sizeof(icer_image_segment_typedef);
                return res;
            }

            data_in_bytes = context.output_ind + (context.output_bit_offset > 0);
            seg->data_length = context.output_ind * 8 + context.output_bit_offset;
            seg->data_crc32 = icer_calculate_segment_crc32(seg);
            seg->crc32 = icer_calculate_packet_crc32(seg);
            output_data->size_used += data_in_bytes;

            segments_encoded[segment_num] = seg;

            segment_num++;
        }
        partition_row_ind += segment_h;
    }

    return ICER_RESULT_OK;
}
#endif

#ifdef USE_DECODE_FUNCTIONS
int icer_decompress_partition_uint8(uint8_t * const data, const partition_param_typdef *params, size_t rowstride,
                                    icer_image_segment_typedef *seg[][7]) {
    int res;
    size_t segment_w, segment_h;
    uint8_t *segment_start;
    uint16_t segment_num = 0;

    size_t partition_col_ind;
    size_t partition_row_ind = 0;

    icer_context_model_typedef context_model;
    icer_decoder_context_typedef context;
    icer_packet_context pkt_context;
    int lsb;
    /*
     * process top region which consists of c columns
     * height of top region is h_t and it contains r_t rows
     */
    for (uint16_t row = 0; row < params->r_t; row++) {
        /*
         * the first r_t0 rows have height y_t
         * the remainder have height y_t + 1
         */
        segment_h = params->y_t + ((row >= params->r_t0) ? 1 : 0);
        partition_col_ind = 0;

        for (uint16_t col = 0; col < params->c; col++) {
            /* the first c_t0 columns have width x_t
             * the remainder have width x_t + 1
             */
            segment_w = params->x_t + ((col >= params->c_t0) ? 1 : 0);
            segment_start = data + partition_row_ind * rowstride + partition_col_ind;
            partition_col_ind += segment_w;

            lsb = ICER_BITPLANES_TO_COMPRESS_8 - 1;
            /* decompress starting from the msb, and stop whenever there is a missing bitplane */
            /* it is impossible to decompress subsequent bit planes if there is a missing bit plane, due to how the context
             * modeller works; it relies on the previously decoded bbitplanes to determine a bit's context */
            while (seg[segment_num][lsb] != NULL && lsb >= 0) {
                pkt_context.subband_type = seg[segment_num][ICER_BITPLANES_TO_COMPRESS_8 - 1]->subband_type;
                pkt_context.lsb = lsb;
                pkt_context.decomp_level = seg[segment_num][ICER_BITPLANES_TO_COMPRESS_8 - 1]->decomp_level;
                icer_init_context_model_vals(&context_model, pkt_context.subband_type);
                icer_init_entropy_decoder_context(&context, (uint8_t *) seg[segment_num][lsb] +
                                                            sizeof(icer_image_segment_typedef),
                                                  seg[segment_num][lsb]->data_length);
                res = icer_decompress_bitplane_uint8(segment_start, segment_w, segment_h, rowstride, &context_model, &context,
                                               &pkt_context);
                if (res != ICER_RESULT_OK) break;
                lsb--;
            }

            segment_num++;
        }
        partition_row_ind += segment_h;
    }

    /*
     * if the bottom region exists, process bottom region
     * which consists of c+1 columns
     */
    for (uint16_t row = 0; row < (params->r - params->r_t); row++) {
        /*
         * the first r_b0 rows have height y_b
         * the remainder have height y_b + 1
         */
        segment_h = params->y_b + ((row >= params->r_b0) ? 1 : 0);
        partition_col_ind = 0;

        for (uint16_t col = 0; col < (params->c + 1); col++) {
            /* the first c_b0 columns have width x_b
             * the remainder have width x_b + 1
             */
            segment_w = params->x_b + ((col >= params->c_b0) ? 1 : 0);
            segment_start = data + partition_row_ind * rowstride + partition_col_ind;
            partition_col_ind += segment_w;

            lsb = ICER_BITPLANES_TO_COMPRESS_8 - 1;
            /* decompress starting from the msb, and stop whenever there is a missing bitplane */
            /* it is impossible to decompress subsequent bit planes if there is a missing bit plane, due to how the context
             * modeller works; it relies on the previously decoded bbitplanes to determine a bit's context */
            while (seg[segment_num][lsb] != NULL && lsb >= 0) {
                pkt_context.subband_type = seg[segment_num][ICER_BITPLANES_TO_COMPRESS_8 - 1]->subband_type;
                pkt_context.lsb = lsb;
                pkt_context.decomp_level = seg[segment_num][ICER_BITPLANES_TO_COMPRESS_8 - 1]->decomp_level;
                icer_init_context_model_vals(&context_model, pkt_context.subband_type);
                icer_init_entropy_decoder_context(&context, (uint8_t *) seg[segment_num][lsb] +
                                                            sizeof(icer_image_segment_typedef),
                                                  seg[segment_num][lsb]->data_length);
                res = icer_decompress_bitplane_uint8(segment_start, segment_w, segment_h, rowstride, &context_model, &context,
                                               &pkt_context);
                if (res != ICER_RESULT_OK) break;
                lsb--;
            }

            segment_num++;
        }
        partition_row_ind += segment_h;
    }

    return ICER_RESULT_OK;
}
#endif
#endif

#ifdef USE_UINT16_FUNCTIONS

#ifdef USE_ENCODE_FUNCTIONS
int icer_compress_partition_uint16(const uint16_t *data, const partition_param_typdef *params, size_t rowstride,
                                   const icer_packet_context *pkt_context, icer_output_data_buf_typedef *output_data,
                                   icer_image_segment_typedef *segments_encoded[]) {
    int res;
    size_t segment_w, segment_h;
    const uint16_t *segment_start;
    uint16_t segment_num = 0;

    size_t partition_col_ind;
    size_t partition_row_ind = 0;

    icer_context_model_typedef context_model;
    icer_encoder_context_typedef context;
    icer_image_segment_typedef *seg;

    uint32_t data_in_bytes;
    /*
     * process top region which consists of c columns
     * height of top region is h_t and it contains r_t rows
     */
    for (uint16_t row = 0; row < params->r_t; row++) {
        /*
         * the first r_t0 rows have height y_t
         * the remainder have height y_t + 1
         */
        segment_h = params->y_t + ((row >= params->r_t0) ? 1 : 0);
        partition_col_ind = 0;

        for (uint16_t col = 0; col < params->c; col++) {
            /* the first c_t0 columns have width x_t
             * the remainder have width x_t + 1
             */
            segment_w = params->x_t + ((col >= params->c_t0) ? 1 : 0);
            segment_start = data + partition_row_ind * rowstride + partition_col_ind;
            partition_col_ind += segment_w;

            icer_init_context_model_vals(&context_model, pkt_context->subband_type);
            res = icer_allocate_data_packet(&seg, output_data, segment_num, pkt_context);
            if (res != ICER_RESULT_OK) return res;

            icer_init_entropy_coder_context(&context, icer_encode_circ_buf, ICER_CIRC_BUF_SIZE,
                                            (uint8_t *) seg + sizeof(icer_image_segment_typedef), seg->data_length);
            res = icer_compress_bitplane_uint16(segment_start, segment_w, segment_h, rowstride, &context_model, &context,
                                               pkt_context);
            if (res != ICER_RESULT_OK) {
                output_data->size_used -= sizeof(icer_image_segment_typedef);
                return res;
            }

            data_in_bytes = context.output_ind + (context.output_bit_offset > 0);
            seg->data_length = context.output_ind * 8 + context.output_bit_offset;
            seg->data_crc32 = icer_calculate_segment_crc32(seg);
            seg->crc32 = icer_calculate_packet_crc32(seg);
            output_data->size_used += data_in_bytes;

            segments_encoded[segment_num] = seg;

            segment_num++;
        }
        partition_row_ind += segment_h;
    }

    /*
     * if the bottom region exists, process bottom region
     * which consists of c+1 columns
     */
    for (uint16_t row = 0; row < (params->r - params->r_t); row++) {
        /*
         * the first r_b0 rows have height y_b
         * the remainder have height y_b + 1
         */
        segment_h = params->y_b + ((row >= params->r_b0) ? 1 : 0);
        partition_col_ind = 0;

        for (uint16_t col = 0; col < (params->c + 1); col++) {
            /* the first c_b0 columns have width x_b
             * the remainder have width x_b + 1
             */
            segment_w = params->x_b + ((col >= params->c_b0) ? 1 : 0);
            segment_start = data + partition_row_ind * rowstride + partition_col_ind;
            partition_col_ind += segment_w;

            icer_init_context_model_vals(&context_model, pkt_context->subband_type);
            res = icer_allocate_data_packet(&seg, output_data, segment_num, pkt_context);
            if (res != ICER_RESULT_OK) return res;

            icer_init_entropy_coder_context(&context, icer_encode_circ_buf, ICER_CIRC_BUF_SIZE,
                                            (uint8_t *) seg + sizeof(icer_image_segment_typedef), seg->data_length);
            res = icer_compress_bitplane_uint16(segment_start, segment_w, segment_h, rowstride, &context_model, &context,
                                               pkt_context);
            if (res != ICER_RESULT_OK) {
                output_data->size_used -= sizeof(icer_image_segment_typedef);
                return res;
            }

            data_in_bytes = context.output_ind + (context.output_bit_offset > 0);
            seg->data_length = context.output_ind * 8 + context.output_bit_offset;
            seg->data_crc32 = icer_calculate_segment_crc32(seg);
            seg->crc32 = icer_calculate_packet_crc32(seg);
            output_data->size_used += data_in_bytes;

            segments_encoded[segment_num] = seg;

            segment_num++;
        }
        partition_row_ind += segment_h;
    }

    return ICER_RESULT_OK;
}
#endif


#ifdef USE_DECODE_FUNCTIONS
int icer_decompress_partition_uint16(uint16_t * const data, const partition_param_typdef * params, size_t rowstride,
                                    icer_image_segment_typedef *seg[][15]) {
    int res;
    size_t segment_w, segment_h;
    uint16_t *segment_start;
    uint16_t segment_num = 0;

    size_t partition_col_ind;
    size_t partition_row_ind = 0;

    icer_context_model_typedef context_model;
    icer_decoder_context_typedef context;
    icer_packet_context pkt_context;
    int lsb;
    /*
     * process top region which consists of c columns
     * height of top region is h_t and it contains r_t rows
     */
    for (uint16_t row = 0; row < params->r_t; row++) {
        /*
         * the first r_t0 rows have height y_t
         * the remainder have height y_t + 1
         */
        segment_h = params->y_t + ((row >= params->r_t0) ? 1 : 0);
        partition_col_ind = 0;

        for (uint16_t col = 0; col < params->c; col++) {
            /* the first c_t0 columns have width x_t
             * the remainder have width x_t + 1
             */
            segment_w = params->x_t + ((col >= params->c_t0) ? 1 : 0);
            segment_start = data + partition_row_ind * rowstride + partition_col_ind;
            partition_col_ind += segment_w;

            lsb = ICER_BITPLANES_TO_COMPRESS_16 - 1;
            /* decompress starting from the msb, and stop whenever there is a missing bitplane */
            /* it is impossible to decompress subsequent bit planes if there is a missing bit plane, due to how the context
             * modeller works; it relies on the previously decoded bbitplanes to determine a bit's context */
            while (seg[segment_num][lsb] != NULL && lsb >= 0) {
                pkt_context.subband_type = seg[segment_num][ICER_BITPLANES_TO_COMPRESS_16 - 1]->subband_type;
                pkt_context.lsb = lsb;
                pkt_context.decomp_level = seg[segment_num][ICER_BITPLANES_TO_COMPRESS_16 - 1]->decomp_level;
                icer_init_context_model_vals(&context_model, pkt_context.subband_type);
                icer_init_entropy_decoder_context(&context, (uint8_t *) seg[segment_num][lsb] +
                                                            sizeof(icer_image_segment_typedef),
                                                  seg[segment_num][lsb]->data_length);
                res = icer_decompress_bitplane_uint16(segment_start, segment_w, segment_h, rowstride, &context_model, &context,
                                                     &pkt_context);
                if (res != ICER_RESULT_OK) break;
                lsb--;
            }

            segment_num++;
        }
        partition_row_ind += segment_h;
    }

    /*
     * if the bottom region exists, process bottom region
     * which consists of c+1 columns
     */
    for (uint16_t row = 0; row < (params->r - params->r_t); row++) {
        /*
         * the first r_b0 rows have height y_b
         * the remainder have height y_b + 1
         */
        segment_h = params->y_b + ((row >= params->r_b0) ? 1 : 0);
        partition_col_ind = 0;

        for (uint16_t col = 0; col < (params->c + 1); col++) {
            /* the first c_b0 columns have width x_b
             * the remainder have width x_b + 1
             */
            segment_w = params->x_b + ((col >= params->c_b0) ? 1 : 0);
            segment_start = data + partition_row_ind * rowstride + partition_col_ind;
            partition_col_ind += segment_w;

            lsb = ICER_BITPLANES_TO_COMPRESS_16 - 1;
            /* decompress starting from the msb, and stop whenever there is a missing bitplane */
            /* it is impossible to decompress subsequent bit planes if there is a missing bit plane, due to how the context
             * modeller works; it relies on the previously decoded bbitplanes to determine a bit's context */
            while (seg[segment_num][lsb] != NULL && lsb >= 0) {
                pkt_context.subband_type = seg[segment_num][ICER_BITPLANES_TO_COMPRESS_16 - 1]->subband_type;
                pkt_context.lsb = lsb;
                pkt_context.decomp_level = seg[segment_num][ICER_BITPLANES_TO_COMPRESS_16 - 1]->decomp_level;
                icer_init_context_model_vals(&context_model, pkt_context.subband_type);
                icer_init_entropy_decoder_context(&context, (uint8_t *) seg[segment_num][lsb] +
                                                            sizeof(icer_image_segment_typedef),
                                                  seg[segment_num][lsb]->data_length);
                res = icer_decompress_bitplane_uint16(segment_start, segment_w, segment_h, rowstride, &context_model, &context,
                                                     &pkt_context);
                if (res != ICER_RESULT_OK) break;
                lsb--;
            }

            segment_num++;
        }
        partition_row_ind += segment_h;
    }

    return ICER_RESULT_OK;
}
#endif
#endif

