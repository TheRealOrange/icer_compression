//
// Created by linyi on 22/3/2023.
//

#include "icer.h"
#include <string.h>

packet_context packets[200];

static inline int comp_packet(const void *a, const void *b) {
    if (((packet_context *)a)->priority == ((packet_context *)b)->priority) {
        return ((packet_context *)a)->subband_type - ((packet_context *)b)->subband_type;
    }
    return ((packet_context *)b)->priority - ((packet_context *)a)->priority;
}

int icer_compress_image_uint8(uint8_t *image, size_t image_w, size_t image_h, uint8_t stages, enum icer_filter_types filt,
                          uint8_t segments, output_data_buf_typedef *output_data) {
    int res;
    res = icer_wavelet_transform_stages(image, image_w, image_h, stages, filt);
    if (res != ICER_RESULT_OK) return res;

    size_t ll_w = get_dim_n_low_stages(image_w, stages);
    size_t ll_h = get_dim_n_low_stages(image_h, stages);

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
        priority = icer_pow_uint(curr_stage, 2);
        for (uint8_t lsb = 0;lsb < 7;lsb++) {
            packets[ind].subband_type = ICER_SUBBAND_HL;
            packets[ind].decomp_level = curr_stage;
            packets[ind].ll_mean_val = ll_mean;
            packets[ind].lsb = lsb;
            packets[ind].priority = priority << lsb;
            packets[ind].image_w = image_w;
            packets[ind].image_h = image_h;
            ind++;

            packets[ind].subband_type = ICER_SUBBAND_LH;
            packets[ind].decomp_level = curr_stage;
            packets[ind].ll_mean_val = ll_mean;
            packets[ind].lsb = lsb;
            packets[ind].priority = priority << lsb;
            packets[ind].image_w = image_w;
            packets[ind].image_h = image_h;
            ind++;

            packets[ind].subband_type = ICER_SUBBAND_HH;
            packets[ind].decomp_level = curr_stage;
            packets[ind].ll_mean_val = ll_mean;
            packets[ind].lsb = lsb;
            packets[ind].priority = (priority/2) << lsb;
            packets[ind].image_w = image_w;
            packets[ind].image_h = image_h;
            ind++;
        }
    }

    priority = icer_pow_uint(stages, 2);
    for (uint8_t lsb = 0;lsb < 7;lsb++) {
        packets[ind].subband_type = ICER_SUBBAND_LL;
        packets[ind].decomp_level = stages;
        packets[ind].ll_mean_val = ll_mean;
        packets[ind].lsb = lsb;
        packets[ind].priority = (2*priority) << lsb;
        packets[ind].image_w = image_w;
        packets[ind].image_h = image_h;
        ind++;
    }

    qsort(packets, ind, sizeof(packet_context), comp_packet);

    partition_param_typdef partition_params;
    uint8_t *data_start;
    for (size_t it = 0;it < ind;it++) {
        if (packets[it].subband_type == ICER_SUBBAND_LL) {
            ll_w = get_dim_n_low_stages(image_w, packets[it].decomp_level);
            ll_h = get_dim_n_low_stages(image_h, packets[it].decomp_level);
            data_start = image;
        } else if (packets[it].subband_type == ICER_SUBBAND_HL) {
            ll_w = get_dim_n_high_stages(image_w, packets[it].decomp_level);
            ll_h = get_dim_n_low_stages(image_h, packets[it].decomp_level);
            data_start = image + get_dim_n_low_stages(image_w, packets[it].decomp_level);
        } else if (packets[it].subband_type == ICER_SUBBAND_LH) {
            ll_w = get_dim_n_low_stages(image_w, packets[it].decomp_level);
            ll_h = get_dim_n_high_stages(image_h, packets[it].decomp_level);
            data_start = image + get_dim_n_low_stages(image_h, packets[it].decomp_level) * image_w;
        } else if (packets[it].subband_type == ICER_SUBBAND_HH) {
            ll_w = get_dim_n_high_stages(image_w, packets[it].decomp_level);
            ll_h = get_dim_n_high_stages(image_h, packets[it].decomp_level);
            data_start = image + get_dim_n_low_stages(image_h, packets[it].decomp_level) * image_w + get_dim_n_low_stages(image_w, packets[it].decomp_level);
        }

        icer_generate_partition_parameters(&partition_params, ll_w, ll_h, segments);
        res = compress_partition_uint8(data_start, &partition_params, image_w, &(packets[it]), output_data);
        if (res != ICER_RESULT_OK) return res;
    }

    return ICER_RESULT_OK;
}

image_segment_typedef *reconstruct_data[ICER_MAX_DECOMP_STAGES+1][ICER_SUBBAND_MAX+1][ICER_MAX_SEGMENTS+1][7];

int icer_decompress_image_uint8(uint8_t *image, size_t *image_w, size_t *image_h, size_t image_bufsize, uint8_t *datastream,
                            size_t data_length, uint8_t stages, enum icer_filter_types filt, uint8_t segments) {
    for (int i = 0;i <= ICER_MAX_DECOMP_STAGES;i++) {
        for (int j = 0;j <= ICER_SUBBAND_MAX;j++) {
            for (int k = 0;k <= ICER_MAX_SEGMENTS;k++) {
                for (int lsb = 0;lsb < 7;lsb++) {
                    reconstruct_data[i][j][k][lsb] = NULL;
                }
            }
        }
    }

    image_segment_typedef *seg = NULL;
    uint8_t *data_start;
    size_t offset = 0;
    size_t res;
    uint16_t ll_mean;
    while ((data_length - offset) > 0) {
        data_start = datastream + offset;
        res = icer_find_packet_in_bytestream(&seg, data_start, data_length - offset);
        if (seg != NULL) {
            offset += res+1;
            reconstruct_data[seg->decomp_level][seg->subband_type][seg->segment_number][seg->lsb] = seg;
            *image_w = seg->image_w;
            *image_h = seg->image_h;
            ll_mean = seg->ll_mean_val;
        } else {
            offset++;
        }
    }

    if (image_bufsize < (*image_w) * (*image_h)) {
        return ICER_BYTE_QUOTA_EXCEEDED;
    }

    size_t ll_w;
    size_t ll_h;
    size_t im_w = *image_w;
    size_t im_h = *image_h;
    memset(image, 0, im_w * im_h * sizeof(image[0]));
    partition_param_typdef partition_params;
    for (uint8_t curr_stage = 1;curr_stage <= stages;curr_stage++) {
        printf("deccomp stage: %d\n", curr_stage);
        if (curr_stage == stages) {
            /* LL subband */
            ll_w = get_dim_n_low_stages(im_w, curr_stage);
            ll_h = get_dim_n_low_stages(im_h, curr_stage);
            data_start = image;

            icer_generate_partition_parameters(&partition_params, ll_w, ll_h, segments);
            decompress_partition_uint8(data_start, &partition_params, im_w, reconstruct_data[curr_stage][ICER_SUBBAND_LL]);
        }

        /* HL subband */
        ll_w = get_dim_n_high_stages(im_w, curr_stage);
        ll_h = get_dim_n_low_stages(im_h, curr_stage);
        data_start = image + get_dim_n_low_stages(im_w, curr_stage);

        icer_generate_partition_parameters(&partition_params, ll_w, ll_h, segments);
        decompress_partition_uint8(data_start, &partition_params, im_w, reconstruct_data[curr_stage][ICER_SUBBAND_HL]);

        /* LH subband */
        ll_w = get_dim_n_low_stages(im_w, curr_stage);
        ll_h = get_dim_n_high_stages(im_h, curr_stage);
        data_start = image + get_dim_n_low_stages(im_h, curr_stage) * im_w;

        icer_generate_partition_parameters(&partition_params, ll_w, ll_h, segments);
        decompress_partition_uint8(data_start, &partition_params, im_w, reconstruct_data[curr_stage][ICER_SUBBAND_LH]);

        /* HH subband */
        ll_w = get_dim_n_high_stages(im_w, curr_stage);
        ll_h = get_dim_n_high_stages(im_h, curr_stage);
        data_start = image + get_dim_n_low_stages(im_h, curr_stage) * im_w + get_dim_n_low_stages(im_w, curr_stage);

        icer_generate_partition_parameters(&partition_params, ll_w, ll_h, segments);
        decompress_partition_uint8(data_start, &partition_params, im_w, reconstruct_data[curr_stage][ICER_SUBBAND_HH]);
    }
    printf("decomp stage done\n");

    icer_from_sign_magnitude_int8(image, im_w * im_h);

    ll_w = get_dim_n_low_stages(im_w, stages);
    ll_h = get_dim_n_low_stages(im_h, stages);
    int8_t *signed_pixel;
    for (size_t row = 0;row < ll_h;row++) {
        signed_pixel = (int8_t*)(image + im_w * row);
        for (size_t col = 0;col < ll_w;col++) {
            (*signed_pixel) += (int8_t)ll_mean;
            signed_pixel++;
        }
    }

    icer_inverse_wavelet_transform_stages(image, im_w, im_h, stages, filt);
    return ICER_RESULT_OK;
}

size_t icer_find_packet_in_bytestream(image_segment_typedef **seg, uint8_t *datastream, size_t data_length) {
    size_t offset = 0;
    (*seg) = NULL;
    while (offset < data_length) {
        (*seg) = (image_segment_typedef*)(datastream + offset);
        if ((*seg)->preamble == ICER_PACKET_PREAMBLE) {
            if ((*seg)->crc32 == icer_calculate_packet_crc32((*seg))) {
               if (icer_ceil_div_uint32((*seg)->data_length, 8) <= (data_length-offset-sizeof(image_segment_typedef))) {
                 return offset;
               }
            } else {
              (*seg) = NULL;
            }
        } else {
          (*seg) = NULL;
        }
        offset++;
    }
    return 0;
}