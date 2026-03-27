#ifndef WIFI_MATH_FORMULA_H
#define WIFI_MATH_FORMULA_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stddef.h>
#include <stdint.h>

typedef struct {
    double mean;
    double sum_sq_diff;
    unsigned long count;
} running_math_stats;

void running_math_stats_init(running_math_stats *rs, uint32_t size);
void update_running_math_stats_param(running_math_stats *rs, double x);
double get_variance_from_math_stats(running_math_stats *rs);
double get_mean_from_math_stats(running_math_stats *rs);
double cal_population_variance(double population_mean, double new_value);

#ifdef __cplusplus
}
#endif

#endif // WIFI_MATH_FORMULA_H
