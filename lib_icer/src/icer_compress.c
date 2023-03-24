//
// Created by linyi on 22/3/2023.
//

#include "icer.h"
#include <string.h>

#ifdef USE_UINT8_FUNCTIONS
icer_packet_context icer_packets[ICER_MAX_PACKETS];
icer_image_segment_typedef *icer_reconstruct_data_8[ICER_MAX_DECOMP_STAGES + 1][ICER_SUBBAND_MAX + 1][ICER_MAX_SEGMENTS + 1][7];
#endif
#ifdef USE_UINT16_FUNCTIONS
icer_packet_context icer_packets_16[ICER_MAX_PACKETS_16];
icer_image_segment_typedef *icer_reconstruct_data_16[ICER_MAX_DECOMP_STAGES + 1][ICER_SUBBAND_MAX + 1][ICER_MAX_SEGMENTS + 1][15];
#endif

static inline int comp_packet(const void *a, const void *b) {
    if (((icer_packet_context *)a)->priority == ((icer_packet_context *)b)->priority) {
        return ((icer_packet_context *)a)->subband_type - ((icer_packet_context *)b)->subband_type;
    }
    return ((icer_packet_context *)b)->priority - ((icer_packet_context *)a)->priority;
}

#ifdef USE_UINT8_FUNCTIONS
int icer_compress_image_uint8(uint8_t *image, size_t image_w, size_t image_h, uint8_t stages, enum icer_filter_types filt,
                              uint8_t segments, icer_output_data_buf_typedef *output_data) {
    int res;
    res = icer_wavelet_transform_stages_uint8(image, image_w, image_h, stages, filt);
    if (res != ICER_RESULT_OK) return res;

    size_t ll_w = icer_get_dim_n_low_stages(image_w, stages);
    size_t ll_h = icer_get_dim_n_low_stages(image_h, stages);

    uint64_t sum = 0;
    uint8_t *pixel;
    for (size_t row = 0;row < ll_h;row++) {
        pixel = image + image_w * row;
        for (size_t col = 0;col < ll_w;col++) {
            sum += (*pixel);
            pixel++;
        }
    }

    uint8_t ll_mean = sum / (ll_w * ll_h);
    if (ll_mean > INT8_MAX) {
        return ICER_INTEGER_OVERFLOW;
    }

    int8_t *signed_pixel;
    for (size_t row = 0;row < ll_h;row++) {
        signed_pixel = (int8_t*)(image + image_w * row);
        for (size_t col = 0;col < ll_w;col++) {
            (*signed_pixel) -= (int8_t)ll_mean;
            signed_pixel++;
        }
    }

    icer_to_sign_magnitude_int8(image, image_w * image_h);

    uint32_t priority = 0;
    uint32_t ind = 0;
    for (uint8_t curr_stage = 1;curr_stage <= stages;curr_stage++) {
        priority = icer_pow_uint(2, curr_stage);
        for (uint8_t lsb = 0;lsb < ICER_BITPLANES_TO_COMPRESS_8;lsb++) {
            icer_packets[ind].subband_type = ICER_SUBBAND_HL;
            icer_packets[ind].decomp_level = curr_stage;
            icer_packets[ind].ll_mean_val = ll_mean;
            icer_packets[ind].lsb = lsb;
            icer_packets[ind].priority = priority << lsb;
            icer_packets[ind].image_w = image_w;
            icer_packets[ind].image_h = image_h;
            ind++; if (ind > ICER_MAX_PACKETS) return ICER_PACKET_COUNT_EXCEEDED;

            icer_packets[ind].subband_type = ICER_SUBBAND_LH;
            icer_packets[ind].decomp_level = curr_stage;
            icer_packets[ind].ll_mean_val = ll_mean;
            icer_packets[ind].lsb = lsb;
            icer_packets[ind].priority = priority << lsb;
            icer_packets[ind].image_w = image_w;
            icer_packets[ind].image_h = image_h;
            ind++; if (ind > ICER_MAX_PACKETS) return ICER_PACKET_COUNT_EXCEEDED;

            icer_packets[ind].subband_type = ICER_SUBBAND_HH;
            icer_packets[ind].decomp_level = curr_stage;
            icer_packets[ind].ll_mean_val = ll_mean;
            icer_packets[ind].lsb = lsb;
            icer_packets[ind].priority = ((priority / 2) << lsb) + 1;
            icer_packets[ind].image_w = image_w;
            icer_packets[ind].image_h = image_h;
            ind++; if (ind > ICER_MAX_PACKETS) return ICER_PACKET_COUNT_EXCEEDED;
        }
    }

    priority = icer_pow_uint(2, stages);
    for (uint8_t lsb = 0;lsb < ICER_BITPLANES_TO_COMPRESS_8;lsb++) {
        icer_packets[ind].subband_type = ICER_SUBBAND_LL;
        icer_packets[ind].decomp_level = stages;
        icer_packets[ind].ll_mean_val = ll_mean;
        icer_packets[ind].lsb = lsb;
        icer_packets[ind].priority = (2 * priority) << lsb;
        icer_packets[ind].image_w = image_w;
        icer_packets[ind].image_h = image_h;
        ind++; if (ind > ICER_MAX_PACKETS) return ICER_PACKET_COUNT_EXCEEDED;

        icer_packets[ind].subband_type = ICER_SUBBAND_LL;
        icer_packets[ind].decomp_level = stages;
        icer_packets[ind].ll_mean_val = ll_mean;
        icer_packets[ind].lsb = lsb;
        icer_packets[ind].priority = priority << lsb;
        icer_packets[ind].image_w = image_w;
        icer_packets[ind].image_h = image_h;
        ind++; if (ind > ICER_MAX_PACKETS) return ICER_PACKET_COUNT_EXCEEDED;
    }

    qsort(icer_packets, ind, sizeof(icer_packet_context), comp_packet);

    partition_param_typdef partition_params;
    uint8_t *data_start = image;
    for (size_t it = 0;it < ind;it++) {
        if (icer_packets[it].subband_type == ICER_SUBBAND_LL) {
            ll_w = icer_get_dim_n_low_stages(image_w, icer_packets[it].decomp_level);
            ll_h = icer_get_dim_n_low_stages(image_h, icer_packets[it].decomp_level);
            data_start = image;
        } else if (icer_packets[it].subband_type == ICER_SUBBAND_HL) {
            ll_w = icer_get_dim_n_high_stages(image_w, icer_packets[it].decomp_level);
            ll_h = icer_get_dim_n_low_stages(image_h, icer_packets[it].decomp_level);
            data_start = image + icer_get_dim_n_low_stages(image_w, icer_packets[it].decomp_level);
        } else if (icer_packets[it].subband_type == ICER_SUBBAND_LH) {
            ll_w = icer_get_dim_n_low_stages(image_w, icer_packets[it].decomp_level);
            ll_h = icer_get_dim_n_high_stages(image_h, icer_packets[it].decomp_level);
            data_start = image + icer_get_dim_n_low_stages(image_h, icer_packets[it].decomp_level) * image_w;
        } else if (icer_packets[it].subband_type == ICER_SUBBAND_HH) {
            ll_w = icer_get_dim_n_high_stages(image_w, icer_packets[it].decomp_level);
            ll_h = icer_get_dim_n_high_stages(image_h, icer_packets[it].decomp_level);
            data_start = image + icer_get_dim_n_low_stages(image_h, icer_packets[it].decomp_level) * image_w +
                         icer_get_dim_n_low_stages(image_w, icer_packets[it].decomp_level);
        }

        icer_generate_partition_parameters(&partition_params, ll_w, ll_h, segments);
        res = icer_compress_partition_uint8(data_start, &partition_params, image_w, &(icer_packets[it]), output_data);
        if (res != ICER_RESULT_OK) {
            return res;
        }
    }

    return ICER_RESULT_OK;
}

int icer_decompress_image_uint8(uint8_t *image, size_t *image_w, size_t *image_h, size_t image_bufsize, uint8_t *datastream,
                            size_t data_length, uint8_t stages, enum icer_filter_types filt, uint8_t segments) {
    for (int i = 0;i <= ICER_MAX_DECOMP_STAGES;i++) {
        for (int j = 0;j <= ICER_SUBBAND_MAX;j++) {
            for (int k = 0;k <= ICER_MAX_SEGMENTS;k++) {
                for (int lsb = 0;lsb < ICER_BITPLANES_TO_COMPRESS_8;lsb++) {
                    icer_reconstruct_data_8[i][j][k][lsb] = NULL;
                }
            }
        }
    }

    icer_image_segment_typedef *seg = NULL;
    uint8_t *seg_start;
    size_t offset = 0;
    size_t res;
    uint16_t ll_mean = 0;
    while ((data_length - offset) > 0) {
        seg_start = datastream + offset;
        res = icer_find_packet_in_bytestream(&seg, seg_start, data_length - offset);
        if (seg != NULL) {
            icer_reconstruct_data_8[seg->decomp_level][seg->subband_type][seg->segment_number][seg->lsb] = seg;
            *image_w = seg->image_w;
            *image_h = seg->image_h;
            ll_mean = seg->ll_mean_val;
        }
        offset += res;
    }

    if (image_bufsize < (*image_w) * (*image_h)) {
        return ICER_BYTE_QUOTA_EXCEEDED;
    }

    uint8_t *data_start;
    size_t ll_w;
    size_t ll_h;
    size_t im_w = *image_w;
    size_t im_h = *image_h;
    memset(image, 0, im_w * im_h * sizeof(image[0]));
    partition_param_typdef partition_params;
    for (uint8_t curr_stage = 1;curr_stage <= stages;curr_stage++) {
        if (curr_stage == stages) {
            /* LL subband */
            ll_w = icer_get_dim_n_low_stages(im_w, curr_stage);
            ll_h = icer_get_dim_n_low_stages(im_h, curr_stage);
            data_start = image;

            res = icer_generate_partition_parameters(&partition_params, ll_w, ll_h, segments);
            if (res != ICER_RESULT_OK) return res;
            res = icer_decompress_partition_uint8(data_start, &partition_params, im_w,
                                                  icer_reconstruct_data_8[curr_stage][ICER_SUBBAND_LL]);
            if (res != ICER_RESULT_OK) return res;
        }

        /* HL subband */
        ll_w = icer_get_dim_n_high_stages(im_w, curr_stage);
        ll_h = icer_get_dim_n_low_stages(im_h, curr_stage);
        data_start = image + icer_get_dim_n_low_stages(im_w, curr_stage);

        res = icer_generate_partition_parameters(&partition_params, ll_w, ll_h, segments);
        if (res != ICER_RESULT_OK) return res;
        res = icer_decompress_partition_uint8(data_start, &partition_params, im_w,
                                              icer_reconstruct_data_8[curr_stage][ICER_SUBBAND_HL]);
        if (res != ICER_RESULT_OK) return res;

        /* LH subband */
        ll_w = icer_get_dim_n_low_stages(im_w, curr_stage);
        ll_h = icer_get_dim_n_high_stages(im_h, curr_stage);
        data_start = image + icer_get_dim_n_low_stages(im_h, curr_stage) * im_w;

        res = icer_generate_partition_parameters(&partition_params, ll_w, ll_h, segments);
        if (res != ICER_RESULT_OK) return res;
        res = icer_decompress_partition_uint8(data_start, &partition_params, im_w,
                                              icer_reconstruct_data_8[curr_stage][ICER_SUBBAND_LH]);
        if (res != ICER_RESULT_OK) return res;

        /* HH subband */
        ll_w = icer_get_dim_n_high_stages(im_w, curr_stage);
        ll_h = icer_get_dim_n_high_stages(im_h, curr_stage);
        data_start = image + icer_get_dim_n_low_stages(im_h, curr_stage) * im_w + icer_get_dim_n_low_stages(im_w, curr_stage);

        res = icer_generate_partition_parameters(&partition_params, ll_w, ll_h, segments);
        if (res != ICER_RESULT_OK) return res;
        res = icer_decompress_partition_uint8(data_start, &partition_params, im_w,
                                              icer_reconstruct_data_8[curr_stage][ICER_SUBBAND_HH]);
        if (res != ICER_RESULT_OK) return res;
    }
    printf("decomp stage done\n");

    icer_from_sign_magnitude_int8(image, im_w * im_h);

    ll_w = icer_get_dim_n_low_stages(im_w, stages);
    ll_h = icer_get_dim_n_low_stages(im_h, stages);
    int8_t *signed_pixel;
    for (size_t row = 0;row < ll_h;row++) {
        signed_pixel = (int8_t*)(image + im_w * row);
        for (size_t col = 0;col < ll_w;col++) {
            (*signed_pixel) += (int8_t)ll_mean;
            signed_pixel++;
        }
    }

    icer_inverse_wavelet_transform_stages_uint8(image, im_w, im_h, stages, filt);
    return ICER_RESULT_OK;
}
#endif

#ifdef USE_UINT8_FUNCTIONS
int icer_compress_image_uint16(uint16_t *image, size_t image_w, size_t image_h, uint8_t stages, enum icer_filter_types filt,
                              uint8_t segments, icer_output_data_buf_typedef *output_data) {
    int res;
    res = icer_wavelet_transform_stages_uint16(image, image_w, image_h, stages, filt);
    if (res != ICER_RESULT_OK) return res;

    size_t ll_w = icer_get_dim_n_low_stages(image_w, stages);
    size_t ll_h = icer_get_dim_n_low_stages(image_h, stages);

    uint64_t sum = 0;
    uint16_t *pixel;
    for (size_t row = 0;row < ll_h;row++) {
        pixel = image + image_w * row;
        for (size_t col = 0;col < ll_w;col++) {
            sum += (*pixel);
            pixel++;
        }
    }

    uint16_t ll_mean = sum / (ll_w * ll_h);
    if (ll_mean > INT16_MAX) {
        return ICER_INTEGER_OVERFLOW;
    }

    int16_t *signed_pixel;
    for (size_t row = 0;row < ll_h;row++) {
        signed_pixel = (int16_t*)(image + image_w * row);
        for (size_t col = 0;col < ll_w;col++) {
            (*signed_pixel) -= (int16_t)ll_mean;
            signed_pixel++;
        }
    }

    icer_to_sign_magnitude_int16(image, image_w * image_h);

    uint64_t priority = 0;
    uint32_t ind = 0;
    for (uint8_t curr_stage = 1;curr_stage <= stages;curr_stage++) {
        priority = icer_pow_uint(2, curr_stage);
        for (uint8_t lsb = 0;lsb < ICER_BITPLANES_TO_COMPRESS_16;lsb++) {
            icer_packets_16[ind].subband_type = ICER_SUBBAND_HL;
            icer_packets_16[ind].decomp_level = curr_stage;
            icer_packets_16[ind].ll_mean_val = ll_mean;
            icer_packets_16[ind].lsb = lsb;
            icer_packets_16[ind].priority = priority << lsb;
            icer_packets_16[ind].image_w = image_w;
            icer_packets_16[ind].image_h = image_h;
            ind++; if (ind > ICER_MAX_PACKETS_16) return ICER_PACKET_COUNT_EXCEEDED;

            icer_packets_16[ind].subband_type = ICER_SUBBAND_LH;
            icer_packets_16[ind].decomp_level = curr_stage;
            icer_packets_16[ind].ll_mean_val = ll_mean;
            icer_packets_16[ind].lsb = lsb;
            icer_packets_16[ind].priority = priority << lsb;
            icer_packets_16[ind].image_w = image_w;
            icer_packets_16[ind].image_h = image_h;
            ind++; if (ind > ICER_MAX_PACKETS_16) return ICER_PACKET_COUNT_EXCEEDED;

            icer_packets_16[ind].subband_type = ICER_SUBBAND_HH;
            icer_packets_16[ind].decomp_level = curr_stage;
            icer_packets_16[ind].ll_mean_val = ll_mean;
            icer_packets_16[ind].lsb = lsb;
            icer_packets_16[ind].priority = ((priority / 2) << lsb) + 1;
            icer_packets_16[ind].image_w = image_w;
            icer_packets_16[ind].image_h = image_h;
            ind++; if (ind > ICER_MAX_PACKETS_16) return ICER_PACKET_COUNT_EXCEEDED;
        }
    }

    priority = icer_pow_uint(2, stages);
    for (uint8_t lsb = 0;lsb < ICER_BITPLANES_TO_COMPRESS_16;lsb++) {
        icer_packets_16[ind].subband_type = ICER_SUBBAND_LL;
        icer_packets_16[ind].decomp_level = stages;
        icer_packets_16[ind].ll_mean_val = ll_mean;
        icer_packets_16[ind].lsb = lsb;
        icer_packets_16[ind].priority = (2 * priority) << lsb;
        icer_packets_16[ind].image_w = image_w;
        icer_packets_16[ind].image_h = image_h;
        ind++; if (ind > ICER_MAX_PACKETS_16) return ICER_PACKET_COUNT_EXCEEDED;

        icer_packets_16[ind].subband_type = ICER_SUBBAND_LL;
        icer_packets_16[ind].decomp_level = stages;
        icer_packets_16[ind].ll_mean_val = ll_mean;
        icer_packets_16[ind].lsb = lsb;
        icer_packets_16[ind].priority = priority << lsb;
        icer_packets_16[ind].image_w = image_w;
        icer_packets_16[ind].image_h = image_h;
        ind++; if (ind > ICER_MAX_PACKETS_16) return ICER_PACKET_COUNT_EXCEEDED;
    }

    qsort(icer_packets_16, ind, sizeof(icer_packet_context), comp_packet);

    partition_param_typdef partition_params;
    uint16_t *data_start = image;
    for (size_t it = 0;it < ind;it++) {
        if (icer_packets_16[it].subband_type == ICER_SUBBAND_LL) {
            ll_w = icer_get_dim_n_low_stages(image_w, icer_packets_16[it].decomp_level);
            ll_h = icer_get_dim_n_low_stages(image_h, icer_packets_16[it].decomp_level);
            data_start = image;
        } else if (icer_packets_16[it].subband_type == ICER_SUBBAND_HL) {
            ll_w = icer_get_dim_n_high_stages(image_w, icer_packets_16[it].decomp_level);
            ll_h = icer_get_dim_n_low_stages(image_h, icer_packets_16[it].decomp_level);
            data_start = image + icer_get_dim_n_low_stages(image_w, icer_packets_16[it].decomp_level);
        } else if (icer_packets_16[it].subband_type == ICER_SUBBAND_LH) {
            ll_w = icer_get_dim_n_low_stages(image_w, icer_packets_16[it].decomp_level);
            ll_h = icer_get_dim_n_high_stages(image_h, icer_packets_16[it].decomp_level);
            data_start = image + icer_get_dim_n_low_stages(image_h, icer_packets_16[it].decomp_level) * image_w;
        } else if (icer_packets_16[it].subband_type == ICER_SUBBAND_HH) {
            ll_w = icer_get_dim_n_high_stages(image_w, icer_packets_16[it].decomp_level);
            ll_h = icer_get_dim_n_high_stages(image_h, icer_packets_16[it].decomp_level);
            data_start = image + icer_get_dim_n_low_stages(image_h, icer_packets_16[it].decomp_level) * image_w +
                         icer_get_dim_n_low_stages(image_w, icer_packets_16[it].decomp_level);
        }

        icer_generate_partition_parameters(&partition_params, ll_w, ll_h, segments);
        res = icer_compress_partition_uint16(data_start, &partition_params, image_w, &(icer_packets_16[it]), output_data);
        if (res != ICER_RESULT_OK) {
            return res;
        }
    }

    return ICER_RESULT_OK;
}

int icer_decompress_image_uint16(uint16_t *image, size_t *image_w, size_t *image_h, size_t image_bufsize, uint8_t *datastream,
                                size_t data_length, uint8_t stages, enum icer_filter_types filt, uint8_t segments) {
    for (int i = 0;i <= ICER_MAX_DECOMP_STAGES;i++) {
        for (int j = 0;j <= ICER_SUBBAND_MAX;j++) {
            for (int k = 0;k <= ICER_MAX_SEGMENTS;k++) {
                for (int lsb = 0;lsb < ICER_BITPLANES_TO_COMPRESS_16;lsb++) {
                    icer_reconstruct_data_16[i][j][k][lsb] = NULL;
                }
            }
        }
    }

    icer_image_segment_typedef *seg = NULL;
    uint8_t *seg_start;
    size_t offset = 0;
    size_t res;
    uint16_t ll_mean = 0;
    while ((data_length - offset) > 0) {
        seg_start = datastream + offset;
        res = icer_find_packet_in_bytestream(&seg, seg_start, data_length - offset);
        if (seg != NULL) {
            icer_reconstruct_data_16[seg->decomp_level][seg->subband_type][seg->segment_number][seg->lsb] = seg;
            *image_w = seg->image_w;
            *image_h = seg->image_h;
            ll_mean = seg->ll_mean_val;
        }
        offset += res;
    }

    if (image_bufsize < (*image_w) * (*image_h)) {
        return ICER_BYTE_QUOTA_EXCEEDED;
    }

    uint16_t *data_start;
    size_t ll_w;
    size_t ll_h;
    size_t im_w = *image_w;
    size_t im_h = *image_h;
    memset(image, 0, im_w * im_h * sizeof(image[0]));
    partition_param_typdef partition_params;
    for (uint8_t curr_stage = 1;curr_stage <= stages;curr_stage++) {
        if (curr_stage == stages) {
            /* LL subband */
            ll_w = icer_get_dim_n_low_stages(im_w, curr_stage);
            ll_h = icer_get_dim_n_low_stages(im_h, curr_stage);
            data_start = image;

            res = icer_generate_partition_parameters(&partition_params, ll_w, ll_h, segments);
            if (res != ICER_RESULT_OK) return res;
            res = icer_decompress_partition_uint16(data_start, &partition_params, im_w,
                                                   icer_reconstruct_data_16[curr_stage][ICER_SUBBAND_LL]);
            if (res != ICER_RESULT_OK) return res;
        }

        /* HL subband */
        ll_w = icer_get_dim_n_high_stages(im_w, curr_stage);
        ll_h = icer_get_dim_n_low_stages(im_h, curr_stage);
        data_start = image + icer_get_dim_n_low_stages(im_w, curr_stage);

        res = icer_generate_partition_parameters(&partition_params, ll_w, ll_h, segments);
        if (res != ICER_RESULT_OK) return res;
        res = icer_decompress_partition_uint16(data_start, &partition_params, im_w,
                                               icer_reconstruct_data_16[curr_stage][ICER_SUBBAND_HL]);
        if (res != ICER_RESULT_OK) return res;

        /* LH subband */
        ll_w = icer_get_dim_n_low_stages(im_w, curr_stage);
        ll_h = icer_get_dim_n_high_stages(im_h, curr_stage);
        data_start = image + icer_get_dim_n_low_stages(im_h, curr_stage) * im_w;

        res = icer_generate_partition_parameters(&partition_params, ll_w, ll_h, segments);
        if (res != ICER_RESULT_OK) return res;
        res = icer_decompress_partition_uint16(data_start, &partition_params, im_w,
                                               icer_reconstruct_data_16[curr_stage][ICER_SUBBAND_LH]);
        if (res != ICER_RESULT_OK) return res;

        /* HH subband */
        ll_w = icer_get_dim_n_high_stages(im_w, curr_stage);
        ll_h = icer_get_dim_n_high_stages(im_h, curr_stage);
        data_start = image + icer_get_dim_n_low_stages(im_h, curr_stage) * im_w + icer_get_dim_n_low_stages(im_w, curr_stage);

        res = icer_generate_partition_parameters(&partition_params, ll_w, ll_h, segments);
        if (res != ICER_RESULT_OK) return res;
        res = icer_decompress_partition_uint16(data_start, &partition_params, im_w,
                                               icer_reconstruct_data_16[curr_stage][ICER_SUBBAND_HH]);
        if (res != ICER_RESULT_OK) return res;
    }
    printf("decomp stage done\n");

    icer_from_sign_magnitude_int16(image, im_w * im_h);

    ll_w = icer_get_dim_n_low_stages(im_w, stages);
    ll_h = icer_get_dim_n_low_stages(im_h, stages);
    int16_t *signed_pixel;
    for (size_t row = 0;row < ll_h;row++) {
        signed_pixel = (int16_t*)(image + im_w * row);
        for (size_t col = 0;col < ll_w;col++) {
            (*signed_pixel) += (int16_t)ll_mean;
            signed_pixel++;
        }
    }

    icer_inverse_wavelet_transform_stages_uint16(image, im_w, im_h, stages, filt);
    return ICER_RESULT_OK;
}
#endif

size_t icer_find_packet_in_bytestream(icer_image_segment_typedef **seg, uint8_t *datastream, size_t data_length) {
    size_t offset = 0;
    (*seg) = NULL;
    while (offset < data_length) {
        (*seg) = (icer_image_segment_typedef*)(datastream + offset);
        if ((*seg)->preamble == ICER_PACKET_PREAMBLE) {
            if ((*seg)->crc32 == icer_calculate_packet_crc32((*seg))) {
               if (icer_ceil_div_uint32((*seg)->data_length, 8) <= (data_length-offset-sizeof(icer_image_segment_typedef))) {
                   if((*seg)->data_crc32 == icer_calculate_segment_crc32((*seg))) {
                       return offset + icer_ceil_div_uint32((*seg)->data_length, 8) + sizeof(icer_image_segment_typedef);
                   }
               }
            }
        }
        (*seg) = NULL;
        offset++;
    }
    return 1;
}