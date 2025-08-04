//
// Created by linyi on 27/3/2023.
//

#include <stdint.h>
#include <string.h>
#include "icer.h"

static inline int comp_packet(const void *a, const void *b) {
    if (((icer_packet_context *)a)->priority == ((icer_packet_context *)b)->priority) {
        return ((icer_packet_context *)a)->subband_type - ((icer_packet_context *)b)->subband_type;
    }
    return ((icer_packet_context *)b)->priority - ((icer_packet_context *)a)->priority;
}

#ifdef USE_UINT8_FUNCTIONS
int icer_compress_image_yuv_uint8(uint8_t *y_channel, uint8_t *u_channel, uint8_t *v_channel, size_t image_w,
                                  size_t image_h, uint8_t stages, enum icer_filter_types filt,
                                  uint8_t segments, icer_output_data_buf_typedef *const output_data) {
    int res;
    res = icer_wavelet_transform_stages_uint8(y_channel, image_w, image_h, stages, filt);
    if (res != ICER_RESULT_OK) return res;

    res = icer_wavelet_transform_stages_uint8(u_channel, image_w, image_h, stages, filt);
    if (res != ICER_RESULT_OK) return res;

    res = icer_wavelet_transform_stages_uint8(v_channel, image_w, image_h, stages, filt);
    if (res != ICER_RESULT_OK) return res;

    size_t ll_w = icer_get_dim_n_low_stages(image_w, stages);
    size_t ll_h = icer_get_dim_n_low_stages(image_h, stages);

    uint64_t sum[ICER_CHANNEL_MAX+1] = { 0 };
    uint8_t *pixel[ICER_CHANNEL_MAX+1];
    for (size_t row = 0;row < ll_h;row++) {
        pixel[ICER_CHANNEL_Y] = y_channel + image_w * row;
        pixel[ICER_CHANNEL_U] = u_channel + image_w * row;
        pixel[ICER_CHANNEL_V] = v_channel + image_w * row;
        for (size_t col = 0;col < ll_w;col++) {
            for (int chan = ICER_CHANNEL_MIN;chan <= ICER_CHANNEL_MAX;chan++) {
                sum[chan] += (*pixel[chan]);
                pixel[chan]++;
            }
        }
    }

    uint8_t ll_mean[ICER_CHANNEL_MAX+1];
    for (int chan = ICER_CHANNEL_MIN;chan <= ICER_CHANNEL_MAX;chan++) {
        ll_mean[chan] = sum[chan] / (ll_w * ll_h);
        if (ll_mean[chan] > INT8_MAX) {
            return ICER_INTEGER_OVERFLOW;
        }
    }

    int8_t *signed_pixel[ICER_CHANNEL_MAX+1];
    for (size_t row = 0;row < ll_h;row++) {
        signed_pixel[ICER_CHANNEL_Y] = (int8_t*)(y_channel + image_w * row);
        signed_pixel[ICER_CHANNEL_U] = (int8_t*)(u_channel + image_w * row);
        signed_pixel[ICER_CHANNEL_V] = (int8_t*)(v_channel + image_w * row);
        for (size_t col = 0;col < ll_w;col++) {
            for (int chan = ICER_CHANNEL_MIN;chan <= ICER_CHANNEL_MAX;chan++) {
                (*signed_pixel[chan]) -= (int8_t)ll_mean[chan];
                signed_pixel[chan] ++;
            }
        }
    }

    icer_to_sign_magnitude_int8(y_channel, image_w * image_h);
    icer_to_sign_magnitude_int8(u_channel, image_w * image_h);
    icer_to_sign_magnitude_int8(v_channel, image_w * image_h);

    uint32_t priority = 0;
    uint32_t ind = 0;
    for (uint8_t curr_stage = 1;curr_stage <= stages;curr_stage++) {
        priority = icer_pow_uint(2, curr_stage);
        for (uint8_t lsb = 0;lsb < ICER_BITPLANES_TO_COMPRESS_8;lsb++) {
            for (int chan = ICER_CHANNEL_MIN;chan <= ICER_CHANNEL_MAX;chan++) {
                if (chan == ICER_CHANNEL_Y) priority *= 2;

                icer_packets[ind].subband_type = ICER_SUBBAND_HL;
                icer_packets[ind].decomp_level = curr_stage;
                icer_packets[ind].ll_mean_val = ll_mean[chan];
                icer_packets[ind].lsb = lsb;
                icer_packets[ind].priority = priority << lsb;
                icer_packets[ind].image_w = image_w;
                icer_packets[ind].image_h = image_h;
                icer_packets[ind].channel = chan;
                ind++; if (ind >= ICER_MAX_PACKETS) return ICER_PACKET_COUNT_EXCEEDED;

                icer_packets[ind].subband_type = ICER_SUBBAND_LH;
                icer_packets[ind].decomp_level = curr_stage;
                icer_packets[ind].ll_mean_val = ll_mean[chan];
                icer_packets[ind].lsb = lsb;
                icer_packets[ind].priority = priority << lsb;
                icer_packets[ind].image_w = image_w;
                icer_packets[ind].image_h = image_h;
                icer_packets[ind].channel = chan;
                ind++; if (ind >= ICER_MAX_PACKETS) return ICER_PACKET_COUNT_EXCEEDED;

                icer_packets[ind].subband_type = ICER_SUBBAND_HH;
                icer_packets[ind].decomp_level = curr_stage;
                icer_packets[ind].ll_mean_val = ll_mean[chan];
                icer_packets[ind].lsb = lsb;
                icer_packets[ind].priority = ((priority / 2) << lsb) + 1;
                icer_packets[ind].image_w = image_w;
                icer_packets[ind].image_h = image_h;
                icer_packets[ind].channel = chan;
                ind++; if (ind >= ICER_MAX_PACKETS) return ICER_PACKET_COUNT_EXCEEDED;
            }
        }
    }

    priority = icer_pow_uint(2, stages);
    for (uint8_t lsb = 0;lsb < ICER_BITPLANES_TO_COMPRESS_8;lsb++) {
        for (int chan = ICER_CHANNEL_MIN;chan <= ICER_CHANNEL_MAX;chan++) {
            if (chan == ICER_CHANNEL_Y) priority *= 2;

            icer_packets[ind].subband_type = ICER_SUBBAND_LL;
            icer_packets[ind].decomp_level = stages;
            icer_packets[ind].ll_mean_val = ll_mean[chan];
            icer_packets[ind].lsb = lsb;
            icer_packets[ind].priority = (2 * priority) << lsb;
            icer_packets[ind].image_w = image_w;
            icer_packets[ind].image_h = image_h;
            icer_packets[ind].channel = chan;
            ind++;
            if (ind >= ICER_MAX_PACKETS) return ICER_PACKET_COUNT_EXCEEDED;
        }

    }

    qsort(icer_packets, ind, sizeof(icer_packet_context), comp_packet);

    for (int i = 0;i <= ICER_MAX_DECOMP_STAGES;i++) {
        for (int j = 0;j <= ICER_SUBBAND_MAX;j++) {
            for (int k = 0;k <= ICER_MAX_SEGMENTS;k++) {
                for (int lsb = 0;lsb < ICER_BITPLANES_TO_COMPRESS_8;lsb++) {
                    for (int chan = ICER_CHANNEL_MIN;chan <= ICER_CHANNEL_MAX;chan++) {
                        icer_rearrange_segments_8[chan][i][j][lsb][k] = NULL;
                    }
                }
            }
        }
    }

    partition_param_typdef partition_params;
    uint8_t *data_chan[ICER_CHANNEL_MAX+1];
    data_chan[ICER_CHANNEL_Y] = y_channel;
    data_chan[ICER_CHANNEL_U] = u_channel;
    data_chan[ICER_CHANNEL_V] = v_channel;
    uint8_t *data_start;
    for (size_t it = 0;it < ind;it++) {
        if (icer_packets[it].subband_type == ICER_SUBBAND_LL) {
            ll_w = icer_get_dim_n_low_stages(image_w, icer_packets[it].decomp_level);
            ll_h = icer_get_dim_n_low_stages(image_h, icer_packets[it].decomp_level);
            data_start = data_chan[icer_packets[it].channel];
        } else if (icer_packets[it].subband_type == ICER_SUBBAND_HL) {
            ll_w = icer_get_dim_n_high_stages(image_w, icer_packets[it].decomp_level);
            ll_h = icer_get_dim_n_low_stages(image_h, icer_packets[it].decomp_level);
            data_start = data_chan[icer_packets[it].channel] + icer_get_dim_n_low_stages(image_w, icer_packets[it].decomp_level);
        } else if (icer_packets[it].subband_type == ICER_SUBBAND_LH) {
            ll_w = icer_get_dim_n_low_stages(image_w, icer_packets[it].decomp_level);
            ll_h = icer_get_dim_n_high_stages(image_h, icer_packets[it].decomp_level);
            data_start = data_chan[icer_packets[it].channel] + icer_get_dim_n_low_stages(image_h, icer_packets[it].decomp_level) * image_w;
        } else if (icer_packets[it].subband_type == ICER_SUBBAND_HH) {
            ll_w = icer_get_dim_n_high_stages(image_w, icer_packets[it].decomp_level);
            ll_h = icer_get_dim_n_high_stages(image_h, icer_packets[it].decomp_level);
            data_start = data_chan[icer_packets[it].channel] + icer_get_dim_n_low_stages(image_h, icer_packets[it].decomp_level) * image_w +
                         icer_get_dim_n_low_stages(image_w, icer_packets[it].decomp_level);
        } else {
            return ICER_FATAL_ERROR;
        }

        icer_generate_partition_parameters(&partition_params, ll_w, ll_h, segments);
        res = icer_compress_partition_uint8(data_start, &partition_params, image_w, &(icer_packets[it]), output_data,
                                            icer_rearrange_segments_8[icer_packets[it].channel][icer_packets[it].decomp_level][icer_packets[it].subband_type][icer_packets[it].lsb]);
        if (res != ICER_RESULT_OK) {
            break;
        }
    }

    size_t rearrange_offset = 0;
    size_t len;
    for (int k = 0;k <= ICER_MAX_SEGMENTS;k++) {
        for (int j = 0;j <= ICER_SUBBAND_MAX;j++) {
            for (int i = 0;i <= ICER_MAX_DECOMP_STAGES;i++) {
                for (int lsb = 0;lsb < ICER_BITPLANES_TO_COMPRESS_8;lsb++) {
                    for (int chan = ICER_CHANNEL_MIN;chan <= ICER_CHANNEL_MAX;chan++) {
                        if (icer_rearrange_segments_8[chan][i][j][lsb][k] != NULL) {
                            len = icer_ceil_div_uint32(icer_rearrange_segments_8[chan][i][j][lsb][k]->data_length, 8) +
                                  sizeof(icer_image_segment_typedef);
                            icer_rearrange_segments_8[chan][i][j][lsb][k]->lsb_chan |= ICER_SET_CHANNEL_MACRO(chan);
                            memcpy(output_data->rearrange_start + rearrange_offset,
                                   icer_rearrange_segments_8[chan][i][j][lsb][k], len);
                            rearrange_offset += len;
                        }
                    }
                }
            }
        }
    }

    return res;
}

int icer_decompress_image_yuv_uint8(uint8_t *y_channel, uint8_t *u_channel, uint8_t *v_channel, size_t *const image_w,
                                    size_t *const image_h, const size_t image_bufsize, const uint8_t *datastream,
                                    const size_t data_length, const uint8_t stages, const enum icer_filter_types filt,
                                    const uint8_t segments) {
    for (int i = 0;i <= ICER_MAX_DECOMP_STAGES;i++) {
        for (int j = 0;j <= ICER_SUBBAND_MAX;j++) {
            for (int k = 0;k <= ICER_MAX_SEGMENTS;k++) {
                for (int lsb = 0;lsb < ICER_BITPLANES_TO_COMPRESS_8;lsb++) {
                    for (int chan = ICER_CHANNEL_MIN;chan <= ICER_CHANNEL_MAX;chan++) {
                        icer_reconstruct_data_8[chan][i][j][k][lsb] = NULL;
                    }
                }
            }
        }
    }

    icer_image_segment_typedef *seg = NULL;
    const uint8_t *seg_start;
    size_t offset = 0;
    size_t pkt_offset;
    int res;
    uint16_t ll_mean[ICER_CHANNEL_MAX + 1];
    while ((data_length - offset) > 0) {
        seg_start = datastream + offset;
        res = icer_find_packet_in_bytestream(&seg, seg_start, data_length - offset, &pkt_offset);
        if (res == ICER_RESULT_OK) {
            icer_reconstruct_data_8[ICER_GET_CHANNEL_MACRO(seg->lsb_chan)][seg->decomp_level][seg->subband_type][seg->segment_number][ICER_GET_LSB_MACRO(seg->lsb_chan)] = seg;
            *image_w = seg->image_w;
            *image_h = seg->image_h;
            ll_mean[ICER_GET_CHANNEL_MACRO(seg->lsb_chan)] = seg->ll_mean_val;
            seg->lsb_chan = ICER_GET_LSB_MACRO(seg->lsb_chan);
        }
        offset += pkt_offset;
    }

    if (image_bufsize < (*image_w) * (*image_h)) {
        return ICER_BYTE_QUOTA_EXCEEDED;
    }

    uint8_t *data_start;
    size_t ll_w;
    size_t ll_h;
    size_t im_w = *image_w;
    size_t im_h = *image_h;
    memset(y_channel, 0, im_w * im_h * sizeof(y_channel[0]));
    memset(u_channel, 0, im_w * im_h * sizeof(u_channel[0]));
    memset(v_channel, 0, im_w * im_h * sizeof(v_channel[0]));
    uint8_t *data_chan[ICER_CHANNEL_MAX + 1];
    data_chan[ICER_CHANNEL_Y] = y_channel;
    data_chan[ICER_CHANNEL_U] = u_channel;
    data_chan[ICER_CHANNEL_V] = v_channel;
    partition_param_typdef partition_params;
    for (uint8_t curr_stage = 1;curr_stage <= stages;curr_stage++) {
        for (int chan = ICER_CHANNEL_MIN;chan <= ICER_CHANNEL_MAX;chan++) {
            if (curr_stage == stages) {
                /* LL subband */
                ll_w = icer_get_dim_n_low_stages(im_w, curr_stage);
                ll_h = icer_get_dim_n_low_stages(im_h, curr_stage);
                data_start = data_chan[chan];

                res = icer_generate_partition_parameters(&partition_params, ll_w, ll_h, segments);
                if (res != ICER_RESULT_OK) return res;
                res = icer_decompress_partition_uint8(data_start, &partition_params, im_w,
                                                      icer_reconstruct_data_8[chan][curr_stage][ICER_SUBBAND_LL]);
                if (res != ICER_RESULT_OK) return res;
            }

            /* HL subband */
            ll_w = icer_get_dim_n_high_stages(im_w, curr_stage);
            ll_h = icer_get_dim_n_low_stages(im_h, curr_stage);
            data_start = data_chan[chan] + icer_get_dim_n_low_stages(im_w, curr_stage);

            res = icer_generate_partition_parameters(&partition_params, ll_w, ll_h, segments);
            if (res != ICER_RESULT_OK) return res;
            res = icer_decompress_partition_uint8(data_start, &partition_params, im_w,
                                                  icer_reconstruct_data_8[chan][curr_stage][ICER_SUBBAND_HL]);
            if (res != ICER_RESULT_OK) return res;

            /* LH subband */
            ll_w = icer_get_dim_n_low_stages(im_w, curr_stage);
            ll_h = icer_get_dim_n_high_stages(im_h, curr_stage);
            data_start = data_chan[chan] + icer_get_dim_n_low_stages(im_h, curr_stage) * im_w;

            res = icer_generate_partition_parameters(&partition_params, ll_w, ll_h, segments);
            if (res != ICER_RESULT_OK) return res;
            res = icer_decompress_partition_uint8(data_start, &partition_params, im_w,
                                                  icer_reconstruct_data_8[chan][curr_stage][ICER_SUBBAND_LH]);
            if (res != ICER_RESULT_OK) return res;

            /* HH subband */
            ll_w = icer_get_dim_n_high_stages(im_w, curr_stage);
            ll_h = icer_get_dim_n_high_stages(im_h, curr_stage);
            data_start = data_chan[chan] + icer_get_dim_n_low_stages(im_h, curr_stage) * im_w +
                         icer_get_dim_n_low_stages(im_w, curr_stage);

            res = icer_generate_partition_parameters(&partition_params, ll_w, ll_h, segments);
            if (res != ICER_RESULT_OK) return res;
            res = icer_decompress_partition_uint8(data_start, &partition_params, im_w,
                                                  icer_reconstruct_data_8[chan][curr_stage][ICER_SUBBAND_HH]);
            if (res != ICER_RESULT_OK) return res;
        }
    }
    printf("decomp stage done\n");

    icer_from_sign_magnitude_int8(y_channel, im_w * im_h);
    icer_from_sign_magnitude_int8(u_channel, im_w * im_h);
    icer_from_sign_magnitude_int8(v_channel, im_w * im_h);

    ll_w = icer_get_dim_n_low_stages(im_w, stages);
    ll_h = icer_get_dim_n_low_stages(im_h, stages);
    int8_t *signed_pixel[ICER_CHANNEL_MAX + 1];
    for (size_t row = 0;row < ll_h;row++) {
        signed_pixel[ICER_CHANNEL_Y] = (int8_t*)(y_channel + im_w * row);
        signed_pixel[ICER_CHANNEL_U] = (int8_t*)(u_channel + im_w * row);
        signed_pixel[ICER_CHANNEL_V] = (int8_t*)(v_channel + im_w * row);
        for (size_t col = 0;col < ll_w;col++) {
            for (int chan = ICER_CHANNEL_MIN;chan <= ICER_CHANNEL_MAX;chan++) {
                (*signed_pixel[chan]) += (int8_t) ll_mean[chan];
                signed_pixel[chan]++;
            }
        }
    }

    icer_inverse_wavelet_transform_stages_uint8(y_channel, im_w, im_h, stages, filt);
    icer_inverse_wavelet_transform_stages_uint8(u_channel, im_w, im_h, stages, filt);
    icer_inverse_wavelet_transform_stages_uint8(v_channel, im_w, im_h, stages, filt);

    icer_remove_negative_uint8(y_channel, im_w, im_h);
    icer_remove_negative_uint8(u_channel, im_w, im_h);
    icer_remove_negative_uint8(v_channel, im_w, im_h);
    return ICER_RESULT_OK;
}
#endif

#ifdef USE_UINT16_FUNCTIONS
#ifdef USE_ENCODE_FUNCTIONS
int icer_compress_image_yuv_uint16(uint16_t *y_channel, uint16_t *u_channel, uint16_t *v_channel, size_t image_w,
                                  size_t image_h, uint8_t stages, enum icer_filter_types filt,
                                  uint8_t segments, icer_output_data_buf_typedef *const output_data) {
    int res;
    res = icer_wavelet_transform_stages_uint16(y_channel, image_w, image_h, stages, filt);
    if (res != ICER_RESULT_OK) return res;

    res = icer_wavelet_transform_stages_uint16(u_channel, image_w, image_h, stages, filt);
    if (res != ICER_RESULT_OK) return res;

    res = icer_wavelet_transform_stages_uint16(v_channel, image_w, image_h, stages, filt);
    if (res != ICER_RESULT_OK) return res;

    size_t ll_w = icer_get_dim_n_low_stages(image_w, stages);
    size_t ll_h = icer_get_dim_n_low_stages(image_h, stages);

    uint64_t sum[ICER_CHANNEL_MAX+1] = { 0 };
    uint16_t *pixel[ICER_CHANNEL_MAX+1];
    for (size_t row = 0;row < ll_h;row++) {
        pixel[ICER_CHANNEL_Y] = y_channel + image_w * row;
        pixel[ICER_CHANNEL_U] = u_channel + image_w * row;
        pixel[ICER_CHANNEL_V] = v_channel + image_w * row;
        for (size_t col = 0;col < ll_w;col++) {
            for (int chan = ICER_CHANNEL_MIN;chan <= ICER_CHANNEL_MAX;chan++) {
                sum[chan] += (*pixel[chan]);
                pixel[chan]++;
            }
        }
    }

    uint16_t ll_mean[ICER_CHANNEL_MAX+1];
    for (int chan = ICER_CHANNEL_MIN;chan <= ICER_CHANNEL_MAX;chan++) {
        ll_mean[chan] = sum[chan] / (ll_w * ll_h);
        if (ll_mean[chan] > INT16_MAX) {
            return ICER_INTEGER_OVERFLOW;
        }
    }

    int16_t *signed_pixel[ICER_CHANNEL_MAX+1];
    for (size_t row = 0;row < ll_h;row++) {
        signed_pixel[ICER_CHANNEL_Y] = (int16_t*)(y_channel + image_w * row);
        signed_pixel[ICER_CHANNEL_U] = (int16_t*)(u_channel + image_w * row);
        signed_pixel[ICER_CHANNEL_V] = (int16_t*)(v_channel + image_w * row);
        for (size_t col = 0;col < ll_w;col++) {
            for (int chan = ICER_CHANNEL_MIN;chan <= ICER_CHANNEL_MAX;chan++) {
                (*signed_pixel[chan]) -= (int16_t)ll_mean[chan];
                signed_pixel[chan] ++;
            }
        }
    }

    icer_to_sign_magnitude_int16(y_channel, image_w * image_h);
    icer_to_sign_magnitude_int16(u_channel, image_w * image_h);
    icer_to_sign_magnitude_int16(v_channel, image_w * image_h);

    uint32_t priority = 0;
    uint32_t ind = 0;
    for (uint8_t curr_stage = 1;curr_stage <= stages;curr_stage++) {
        priority = icer_pow_uint(2, curr_stage);
        for (uint8_t lsb = 0;lsb < ICER_BITPLANES_TO_COMPRESS_16;lsb++) {
            for (int chan = ICER_CHANNEL_MIN;chan <= ICER_CHANNEL_MAX;chan++) {
                if (chan == ICER_CHANNEL_Y) priority *= 2;

                icer_packets_16[ind].subband_type = ICER_SUBBAND_HL;
                icer_packets_16[ind].decomp_level = curr_stage;
                icer_packets_16[ind].ll_mean_val = ll_mean[chan];
                icer_packets_16[ind].lsb = lsb;
                icer_packets_16[ind].priority = priority << lsb;
                icer_packets_16[ind].image_w = image_w;
                icer_packets_16[ind].image_h = image_h;
                icer_packets_16[ind].channel = chan;
                ind++; if (ind >= ICER_MAX_PACKETS_16) return ICER_PACKET_COUNT_EXCEEDED;

                icer_packets_16[ind].subband_type = ICER_SUBBAND_LH;
                icer_packets_16[ind].decomp_level = curr_stage;
                icer_packets_16[ind].ll_mean_val = ll_mean[chan];
                icer_packets_16[ind].lsb = lsb;
                icer_packets_16[ind].priority = priority << lsb;
                icer_packets_16[ind].image_w = image_w;
                icer_packets_16[ind].image_h = image_h;
                icer_packets_16[ind].channel = chan;
                ind++; if (ind >= ICER_MAX_PACKETS_16) return ICER_PACKET_COUNT_EXCEEDED;

                icer_packets_16[ind].subband_type = ICER_SUBBAND_HH;
                icer_packets_16[ind].decomp_level = curr_stage;
                icer_packets_16[ind].ll_mean_val = ll_mean[chan];
                icer_packets_16[ind].lsb = lsb;
                icer_packets_16[ind].priority = ((priority / 2) << lsb) + 1;
                icer_packets_16[ind].image_w = image_w;
                icer_packets_16[ind].image_h = image_h;
                icer_packets_16[ind].channel = chan;
                ind++; if (ind >= ICER_MAX_PACKETS_16) return ICER_PACKET_COUNT_EXCEEDED;
            }
        }
    }

    priority = icer_pow_uint(2, stages);
    for (uint8_t lsb = 0;lsb < ICER_BITPLANES_TO_COMPRESS_16;lsb++) {
        for (int chan = ICER_CHANNEL_MIN;chan <= ICER_CHANNEL_MAX;chan++) {
            if (chan == ICER_CHANNEL_Y) priority *= 2;

            icer_packets_16[ind].subband_type = ICER_SUBBAND_LL;
            icer_packets_16[ind].decomp_level = stages;
            icer_packets_16[ind].ll_mean_val = ll_mean[chan];
            icer_packets_16[ind].lsb = lsb;
            icer_packets_16[ind].priority = (2 * priority) << lsb;
            icer_packets_16[ind].image_w = image_w;
            icer_packets_16[ind].image_h = image_h;
            icer_packets_16[ind].channel = chan;
            ind++;
            if (ind >= ICER_MAX_PACKETS_16) return ICER_PACKET_COUNT_EXCEEDED;
        }

    }

    qsort(icer_packets_16, ind, sizeof(icer_packet_context), comp_packet);

    for (int i = 0;i <= ICER_MAX_DECOMP_STAGES;i++) {
        for (int j = 0;j <= ICER_SUBBAND_MAX;j++) {
            for (int k = 0;k <= ICER_MAX_SEGMENTS;k++) {
                for (int lsb = 0;lsb < ICER_BITPLANES_TO_COMPRESS_16;lsb++) {
                    for (int chan = ICER_CHANNEL_MIN;chan <= ICER_CHANNEL_MAX;chan++) {
                        icer_rearrange_segments_16[chan][i][j][lsb][k] = NULL;
                    }
                }
            }
        }
    }

    partition_param_typdef partition_params;
    uint16_t *data_chan[ICER_CHANNEL_MAX+1];
    data_chan[ICER_CHANNEL_Y] = y_channel;
    data_chan[ICER_CHANNEL_U] = u_channel;
    data_chan[ICER_CHANNEL_V] = v_channel;
    uint16_t *data_start;
    for (size_t it = 0;it < ind;it++) {
        if (icer_packets_16[it].subband_type == ICER_SUBBAND_LL) {
            ll_w = icer_get_dim_n_low_stages(image_w, icer_packets_16[it].decomp_level);
            ll_h = icer_get_dim_n_low_stages(image_h, icer_packets_16[it].decomp_level);
            data_start = data_chan[icer_packets_16[it].channel];
        } else if (icer_packets_16[it].subband_type == ICER_SUBBAND_HL) {
            ll_w = icer_get_dim_n_high_stages(image_w, icer_packets_16[it].decomp_level);
            ll_h = icer_get_dim_n_low_stages(image_h, icer_packets_16[it].decomp_level);
            data_start = data_chan[icer_packets_16[it].channel] + icer_get_dim_n_low_stages(image_w, icer_packets_16[it].decomp_level);
        } else if (icer_packets_16[it].subband_type == ICER_SUBBAND_LH) {
            ll_w = icer_get_dim_n_low_stages(image_w, icer_packets_16[it].decomp_level);
            ll_h = icer_get_dim_n_high_stages(image_h, icer_packets_16[it].decomp_level);
            data_start = data_chan[icer_packets_16[it].channel] + icer_get_dim_n_low_stages(image_h, icer_packets_16[it].decomp_level) * image_w;
        } else if (icer_packets_16[it].subband_type == ICER_SUBBAND_HH) {
            ll_w = icer_get_dim_n_high_stages(image_w, icer_packets_16[it].decomp_level);
            ll_h = icer_get_dim_n_high_stages(image_h, icer_packets_16[it].decomp_level);
            data_start = data_chan[icer_packets_16[it].channel] + icer_get_dim_n_low_stages(image_h, icer_packets_16[it].decomp_level) * image_w +
                         icer_get_dim_n_low_stages(image_w, icer_packets_16[it].decomp_level);
        } else {
            return ICER_FATAL_ERROR;
        }

        icer_generate_partition_parameters(&partition_params, ll_w, ll_h, segments);
        res = icer_compress_partition_uint16(data_start, &partition_params, image_w, &(icer_packets_16[it]), output_data,
                                            icer_rearrange_segments_16[icer_packets_16[it].channel][icer_packets_16[it].decomp_level][icer_packets_16[it].subband_type][icer_packets_16[it].lsb]);
        if (res != ICER_RESULT_OK) {
            break;
        }
    }

    size_t rearrange_offset = 0;
    size_t len;
    for (int k = 0;k <= ICER_MAX_SEGMENTS;k++) {
        for (int j = ICER_SUBBAND_MAX;j >= 0;j--) {
            for (int i = ICER_MAX_DECOMP_STAGES;i >= 0;i--) {
                for (int lsb = ICER_BITPLANES_TO_COMPRESS_16 - 1;lsb >= 0;lsb--) {
                    for (int chan = ICER_CHANNEL_MIN;chan <= ICER_CHANNEL_MAX;chan++) {
                        if (icer_rearrange_segments_16[chan][i][j][lsb][k] != NULL) {
                            len = icer_ceil_div_uint32(icer_rearrange_segments_16[chan][i][j][lsb][k]->data_length, 8) +
                                  sizeof(icer_image_segment_typedef);
                            icer_rearrange_segments_16[chan][i][j][lsb][k]->lsb_chan |= ICER_SET_CHANNEL_MACRO(chan);
                            memcpy(output_data->rearrange_start + rearrange_offset,
                                   icer_rearrange_segments_16[chan][i][j][lsb][k], len);
                            rearrange_offset += len;
                        }
                    }
                }
            }
        }
    }

    return res;
}
#endif

#ifdef USE_DECODE_FUNCTIONS
int icer_decompress_image_yuv_uint16(uint16_t * const y_channel, uint16_t * const u_channel, uint16_t * const v_channel, size_t *const image_w,
                                    size_t *const image_h, const size_t image_bufsize, const uint8_t *datastream,
                                    const size_t data_length, const uint8_t stages, const enum icer_filter_types filt,
                                    const uint8_t segments) {
    for (int i = 0;i <= ICER_MAX_DECOMP_STAGES;i++) {
        for (int j = 0;j <= ICER_SUBBAND_MAX;j++) {
            for (int k = 0;k <= ICER_MAX_SEGMENTS;k++) {
                for (int lsb = 0;lsb < ICER_BITPLANES_TO_COMPRESS_16;lsb++) {
                    for (int chan = ICER_CHANNEL_MIN;chan <= ICER_CHANNEL_MAX;chan++) {
                        icer_reconstruct_data_16[chan][i][j][k][lsb] = NULL;
                    }
                }
            }
        }
    }

    icer_image_segment_typedef *seg = NULL;
    const uint8_t *seg_start;
    size_t offset = 0;
    size_t pkt_offset;
    int res;
    uint16_t ll_mean[ICER_CHANNEL_MAX + 1];
    while ((data_length - offset) > 0) {
        seg_start = datastream + offset;
        res = icer_find_packet_in_bytestream(&seg, seg_start, data_length - offset, &pkt_offset);
        if (res == ICER_RESULT_OK) {
            icer_reconstruct_data_16[ICER_GET_CHANNEL_MACRO(seg->lsb_chan)][seg->decomp_level][seg->subband_type][seg->segment_number][ICER_GET_LSB_MACRO(seg->lsb_chan)] = seg;
            *image_w = seg->image_w;
            *image_h = seg->image_h;
            ll_mean[ICER_GET_CHANNEL_MACRO(seg->lsb_chan)] = seg->ll_mean_val;
            seg->lsb_chan = ICER_GET_LSB_MACRO(seg->lsb_chan);
        }
        offset += pkt_offset;
    }

    if (image_bufsize < (*image_w) * (*image_h)) {
        return ICER_BYTE_QUOTA_EXCEEDED;
    }

    uint16_t *data_start;
    size_t ll_w;
    size_t ll_h;
    size_t im_w = *image_w;
    size_t im_h = *image_h;
    memset(y_channel, 0, im_w * im_h * sizeof(y_channel[0]));
    memset(u_channel, 0, im_w * im_h * sizeof(u_channel[0]));
    memset(v_channel, 0, im_w * im_h * sizeof(v_channel[0]));
    uint16_t *data_chan[ICER_CHANNEL_MAX + 1];
    data_chan[ICER_CHANNEL_Y] = y_channel;
    data_chan[ICER_CHANNEL_U] = u_channel;
    data_chan[ICER_CHANNEL_V] = v_channel;
    partition_param_typdef partition_params;
    for (uint8_t curr_stage = 1;curr_stage <= stages;curr_stage++) {
        for (int chan = ICER_CHANNEL_MIN;chan <= ICER_CHANNEL_MAX;chan++) {
            if (curr_stage == stages) {
                /* LL subband */
                ll_w = icer_get_dim_n_low_stages(im_w, curr_stage);
                ll_h = icer_get_dim_n_low_stages(im_h, curr_stage);
                data_start = data_chan[chan];

                res = icer_generate_partition_parameters(&partition_params, ll_w, ll_h, segments);
                if (res != ICER_RESULT_OK) return res;
                res = icer_decompress_partition_uint16(data_start, &partition_params, im_w,
                                                      icer_reconstruct_data_16[chan][curr_stage][ICER_SUBBAND_LL]);
                if (res != ICER_RESULT_OK) return res;
            }

            /* HL subband */
            ll_w = icer_get_dim_n_high_stages(im_w, curr_stage);
            ll_h = icer_get_dim_n_low_stages(im_h, curr_stage);
            data_start = data_chan[chan] + icer_get_dim_n_low_stages(im_w, curr_stage);

            res = icer_generate_partition_parameters(&partition_params, ll_w, ll_h, segments);
            if (res != ICER_RESULT_OK) return res;
            res = icer_decompress_partition_uint16(data_start, &partition_params, im_w,
                                                  icer_reconstruct_data_16[chan][curr_stage][ICER_SUBBAND_HL]);
            if (res != ICER_RESULT_OK) return res;

            /* LH subband */
            ll_w = icer_get_dim_n_low_stages(im_w, curr_stage);
            ll_h = icer_get_dim_n_high_stages(im_h, curr_stage);
            data_start = data_chan[chan] + icer_get_dim_n_low_stages(im_h, curr_stage) * im_w;

            res = icer_generate_partition_parameters(&partition_params, ll_w, ll_h, segments);
            if (res != ICER_RESULT_OK) return res;
            res = icer_decompress_partition_uint16(data_start, &partition_params, im_w,
                                                  icer_reconstruct_data_16[chan][curr_stage][ICER_SUBBAND_LH]);
            if (res != ICER_RESULT_OK) return res;

            /* HH subband */
            ll_w = icer_get_dim_n_high_stages(im_w, curr_stage);
            ll_h = icer_get_dim_n_high_stages(im_h, curr_stage);
            data_start = data_chan[chan] + icer_get_dim_n_low_stages(im_h, curr_stage) * im_w +
                         icer_get_dim_n_low_stages(im_w, curr_stage);

            res = icer_generate_partition_parameters(&partition_params, ll_w, ll_h, segments);
            if (res != ICER_RESULT_OK) return res;
            res = icer_decompress_partition_uint16(data_start, &partition_params, im_w,
                                                  icer_reconstruct_data_16[chan][curr_stage][ICER_SUBBAND_HH]);
            if (res != ICER_RESULT_OK) return res;
        }
    }

    icer_from_sign_magnitude_int16(y_channel, im_w * im_h);
    icer_from_sign_magnitude_int16(u_channel, im_w * im_h);
    icer_from_sign_magnitude_int16(v_channel, im_w * im_h);

    ll_w = icer_get_dim_n_low_stages(im_w, stages);
    ll_h = icer_get_dim_n_low_stages(im_h, stages);
    int16_t *signed_pixel[ICER_CHANNEL_MAX + 1];
    for (size_t row = 0;row < ll_h;row++) {
        signed_pixel[ICER_CHANNEL_Y] = (int16_t*)(y_channel + im_w * row);
        signed_pixel[ICER_CHANNEL_U] = (int16_t*)(u_channel + im_w * row);
        signed_pixel[ICER_CHANNEL_V] = (int16_t*)(v_channel + im_w * row);
        for (size_t col = 0;col < ll_w;col++) {
            for (int chan = ICER_CHANNEL_MIN;chan <= ICER_CHANNEL_MAX;chan++) {
                (*signed_pixel[chan]) += (int16_t) ll_mean[chan];
                signed_pixel[chan]++;
            }
        }
    }

    icer_inverse_wavelet_transform_stages_uint16(y_channel, im_w, im_h, stages, filt);
    icer_inverse_wavelet_transform_stages_uint16(u_channel, im_w, im_h, stages, filt);
    icer_inverse_wavelet_transform_stages_uint16(v_channel, im_w, im_h, stages, filt);

    icer_remove_negative_uint16(y_channel, im_w, im_h);
    icer_remove_negative_uint16(u_channel, im_w, im_h);
    icer_remove_negative_uint16(v_channel, im_w, im_h);
    return ICER_RESULT_OK;
}
#endif
#endif