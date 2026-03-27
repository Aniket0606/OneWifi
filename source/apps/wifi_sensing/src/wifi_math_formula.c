#include "wifi_math_formula.h"
#include <stdio.h>
#include <string.h>

void running_math_stats_init(running_math_stats *rs, uint32_t size)
{
    for (uint32_t index = 0; index < size; index++) {
        rs[index].mean = 0.0;
        rs[index].sum_sq_diff = 0.0;
        rs[index].count = 0;
    }
}

void update_running_math_stats_param(running_math_stats *rs, double x)
{
    rs->count++;

    double delta = x - rs->mean;
    rs->mean += delta / rs->count;
    double delta2 = x - rs->mean;

    rs->sum_sq_diff += delta * delta2;
}

double get_variance_from_math_stats(running_math_stats *rs)
{
    if (rs->count == 0) {
        return 0.0;
    }
    return rs->sum_sq_diff / rs->count;   // population variance
}

double get_mean_from_math_stats(running_math_stats *rs)
{
    return rs->mean;
}

double cal_population_variance(double population_mean, double new_value)
{
    return ((new_value - population_mean) / 2);
}
