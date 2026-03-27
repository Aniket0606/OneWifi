#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include "wifi_hal.h"
#include "wifi_util.h"
#include "wifi_math_formula.h"
#include "motion_core.h"

void set_cal_duration_for_each_clients(uint32_t duration_sec)
{
    motion_core_param_t *p_cfg = (motion_core_param_t *)get_motion_core_param();
    uint32_t cal_packets = ((duration_sec * 1000) / CSI_MOTION_CORE_INTERVAL);
    motion_whitelist_info_t *p_sta_info;

    p_sta_info = (motion_whitelist_info_t *)hash_map_get_first(p_cfg->motion_sta_map);
    while (p_sta_info != NULL) {
        //we need to add lock for this variable protection
        p_sta_info->cal_packets_cnt = cal_packets;
        wifi_util_info_print(WIFI_MGR,"%s:%d number of packet:%d calibration"
            "configure for sta:%s\r\n", __func__, __LINE__,
            cal_packets, p_sta_info->sta_mac);
        p_sta_info = hash_map_get_next(p_cfg->motion_sta_map, p_sta_info);
    }
}

static inline double calculate_magnitude(double l_real, double l_imaginary)
{
    return sqrt(l_real * l_real + l_imaginary * l_imaginary);
}

static inline void convert_re_and_im(uint32_t per_tone_data_32bit, double *real, double *img)
{
    uint16_t real_u = per_tone_data_32bit & 0xFFFF;
    uint16_t imag_u = (per_tone_data_32bit >> 16) & 0xFFFF;

    // Convert unsigned to signed 16-bit (two's complement)
    int16_t real_s = (real_u >= 0x8000) ? (int16_t)(real_u - 0x10000) : (int16_t)real_u;
    int16_t imag_s = (imag_u >= 0x8000) ? (int16_t)(imag_u - 0x10000) : (int16_t)imag_u;

    // Convert S9.6 fixed-point to float by dividing by 64
    *real = real_s / 64.0f;
    *img = imag_s / 64.0f;
}

int process_csi_motion_data(char *str_mac_addr, wifi_csi_data_t  *p_csi,
    motion_core_param_t *p_core_cfg, motion_whitelist_info_t *p_cfg)
{
    double cal_variance, cal_mean;
    double avg_magnitude = 0.0, magnitude = 0.0;
    uint32_t ant_index, sub_car_index;
    double real, img;

    for (ant_index = 0; ant_index < p_csi->frame_info.Nr; ant_index++) {
        avg_magnitude = 0.0;
        for (sub_car_index = 0; sub_car_index < p_csi->frame_info.num_sc; sub_car_index++) {
            convert_re_and_im(p_csi->csi_matrix[sub_car_index][ant_index][0], &real, &img);
            magnitude = calculate_magnitude(real, img);
            avg_magnitude += magnitude;
        }
        avg_magnitude /= p_csi->frame_info.num_sc;

        if (p_cfg->cal_packets_cnt > 0) {
            update_running_math_stats_param(&p_cfg->csi_cal_stats_param[ant_index], avg_magnitude);
            cal_variance = get_variance_from_math_stats(&p_cfg->csi_cal_stats_param[ant_index]);
            cal_mean = get_mean_from_math_stats(&p_cfg->csi_cal_stats_param[ant_index]);
            wifi_util_dbg_print(WIFI_MGR,"===>Ant:%d avg_magnitude:%f variance:%f mean:%f\r\n",
                ant_index, avg_magnitude, cal_variance, cal_mean);
        } else {
            cal_variance = get_variance_from_math_stats(&p_cfg->csi_cal_stats_param[ant_index]);
            cal_mean = get_mean_from_math_stats(&p_cfg->csi_cal_stats_param[ant_index]);
            double l_variance = cal_population_variance(cal_mean, avg_magnitude);
            wifi_util_info_print(WIFI_MGR,"===>Ant:%d avg_magnitude:%f variance:%f (diff:%f)\r\n",
                ant_index, avg_magnitude, l_variance, (cal_variance - l_variance));
        }
    }

    if (p_cfg->cal_packets_cnt > 0) {
        p_cfg->cal_packets_cnt--;
    }

    return 0;
}
