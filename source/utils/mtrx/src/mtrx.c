/**
 * Copyright 2023 Comcast Cable Communications Management, LLC
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdbool.h>
#include <errno.h>
#include <signal.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>
#include "mtrx.h"
#include "wifi_util.h"

int multiply_polynomials(vector_t *out, vector_t *v1, vector_t *v2)
{
    matrix_t m1, m2, m3;
    unsigned int i, j;

    m1.rows = v1->num + v2->num - 1;
    m1.cols = v2->num;

    nullify(&m1);

    for (j = 0; j < m1.cols; j++) {
        for (i = j; i < m1.rows; i++) {
            m1.val[i][j] = v1->val[i - j];
        }
    }

    m2.rows = v2->num;
    m2.cols = 1;

    for (i = 0; i < m2.rows; i++) {
        m2.val[i][0] = v2->val[i];
    }

    print_matrix(&m1);
    print_matrix(&m2);

    if (multiply(&m3, &m1, &m2) != 0) {
        wifi_util_error_print(WIFI_CSI,"%s:%d: Multiplcation can not be done\n", __func__, __LINE__);
        return -1;
    }
    print_matrix(&m3);

    out->num = m3.rows;

    for (i = 0; i < m3.rows; i++) {
        out->val[i] = m3.val[i][0];
    }

    return 0;
}

int nullify_s(matrix_s_t *in)
{
    unsigned int i, j;

    for (i = 0; i < in->rows; i++) {
        for (j = 0; j < in->cols; j++) {
            memset(in->val[i][j], 0, sizeof(expression_t));
        }
    }
    
    return 0;
}

int nullify(matrix_t *in)
{
    unsigned int i, j;

    for (i = 0; i < in->rows; i++) {
        for (j = 0; j < in->cols; j++) {
            in->val[i][j].re = 0.0;
        }
    }

    return 0;
}

int adjoint(matrix_t *out, matrix_t *in)
{
    matrix_t t, m;
    double det = 0.0;
    unsigned int i, j;

    if (in->rows != in->cols) {
        wifi_util_error_print(WIFI_CSI,"%s:%d: Can not find inverse of matrix\n", __func__, __LINE__);
        return -1;
    }

    if (transpose(&t, in) != 0) {
        wifi_util_error_print(WIFI_CSI,"%s:%d: Can not find transpose\n", __func__, __LINE__);
        return -1;
    }

    for (i = 0; i < t.rows; i++) {
        for (j = 0; j < t.cols; j++) {
            minor(&m, &t, i, j);
            if (determinant(&det, &m) != 0) {
                wifi_util_error_print(WIFI_CSI,"%s:%d: Can not find determinant of matrix\n", __func__, __LINE__);
                return -1;
            }

            out->val[i][j].re = pow(-1, i + j) * det;
            nullify(&m);
        }
    }

    out->rows = in->rows;
    out->cols = in->cols;

    return 0;
}

int inverse(matrix_t *out, matrix_t *in)
{
    double det;
    matrix_t a;
    unsigned int i, j;

    if (in->rows != in->cols) {
        wifi_util_error_print(WIFI_CSI,"%s:%d: Can not find inverse of matrix\n", __func__, __LINE__);
        return -1;
    }

    if (determinant(&det, in) != 0) {
        wifi_util_error_print(WIFI_CSI,"%s:%d: Can not find inverse of matrix\n", __func__, __LINE__);
        return -1;
    }

    if (det == 0) {
        wifi_util_error_print(WIFI_CSI,"%s:%d: Matrix is singular\n", __func__, __LINE__);
        return -1;
    }

    if (adjoint(&a, in) != 0) {
        wifi_util_error_print(WIFI_CSI,"%s:%d: Cannot find adjoint\n", __func__, __LINE__);
        return -1;
    }

    out->rows = in->rows;
    out->cols = in->cols;

    for (i = 0; i < in->rows; i++) {
        for (j = 0; j < in->cols; j++) {
            out->val[i][j].re = a.val[i][j].re/det;
        }
    }

    return 0;
}

int clone(matrix_t *out, matrix_t *in)
{
    unsigned int i, j;

    out->rows = in->rows;
    out->cols = in->cols;

    for (i = 0; i < in->rows; i++) {
        for (j = 0; j < in->cols; j++) {
            out->val[i][j] = in->val[i][j];
        }
    }

    return 0;
}

int minor_s(matrix_s_t *out, matrix_s_t *in, unsigned int row, unsigned int col)
{
    unsigned int i, j, p = 0, q = 0;

    out->rows = in->rows - 1;
    out->cols = in->cols - 1;

    nullify_s(out);

    for (i = 0; i < in->rows; i++) {
        if (i == row) {
            continue;
        }
        for (j = 0; j < in->cols; j++) {
            if (j == col) {
                continue;
            }
            strncpy(out->val[p][q], in->val[i][j], strlen(in->val[i][j]) + 1); q++;
        }
        q = 0;
        p++;
    }

    return 0;
}

int minor(matrix_t *out, matrix_t *in, unsigned int row, unsigned int col)
{
    unsigned int i, j, p = 0, q = 0;

    out->rows = in->rows - 1;
    out->cols = in->cols - 1;

    for (i = 0; i < in->rows; i++) {
        if (i == row) {
            continue;
        }
        for (j = 0; j < in->cols; j++) {
            if (j == col) {
                continue;
            }
            out->val[p][q] = in->val[i][j]; q++;
        }
        q = 0;
        p++;
    }

    return 0;
}

int arguments_vector(vector_t *out, large_expression_t e)
{
    large_expression_t copy = {0};
    char *tmp, *s;
    unsigned int i, degree, highest = 0;
    vector_s_t  v_s_tmp = {0}, v_s_signed = {0};
    
    strncpy(copy, e, strlen(e) + 1);
    tmp = copy;
    s = copy;
    
    v_s_tmp.num = 0;
    v_s_signed.num = 0;
    
    if (*tmp == '-') {
        strncpy(v_s_tmp.val[v_s_tmp.num], "-", strlen("-") + 1);
        v_s_tmp.num++;
        tmp++;
    } else if (*tmp == '+') {
        strncpy(v_s_tmp.val[v_s_tmp.num], "+", strlen("+") + 1);
        v_s_tmp.num++;
        tmp++;
    } else {
        strncpy(v_s_tmp.val[v_s_tmp.num], "+", strlen("+") + 1);
        v_s_tmp.num++;
    }
            
    while (tmp != NULL) {
        if (*tmp == '+') {
            strncpy(v_s_tmp.val[v_s_tmp.num], "+", strlen("+") + 1);
            v_s_tmp.num++;
            
        } else if (*tmp == '-') {
            strncpy(v_s_tmp.val[v_s_tmp.num], "-", strlen("-") + 1);
            v_s_tmp.num++;
            
        } else if (*tmp != ' ') {
            s = tmp;
            if ((tmp = strchr(s, ' ')) != NULL) {
                *tmp = 0;
                strncpy(v_s_tmp.val[v_s_tmp.num], s, strlen(s) + 1);
                v_s_tmp.num++;
            } else if ((tmp = strchr(s, '+')) != NULL) {
                *tmp = 0;
                strncpy(v_s_tmp.val[v_s_tmp.num], s, strlen(s) + 1);
                v_s_tmp.num++;
                strncpy(v_s_tmp.val[v_s_tmp.num], "+", strlen("+") + 1);
                v_s_tmp.num++;
            } else if ((tmp = strchr(s, '-')) != NULL) {
                *tmp = 0;
                strncpy(v_s_tmp.val[v_s_tmp.num], s, strlen(s) + 1);
                v_s_tmp.num++;
                strncpy(v_s_tmp.val[v_s_tmp.num], "-", strlen("-") + 1);
                v_s_tmp.num++;
            } else if (tmp == NULL) {
                strncpy(v_s_tmp.val[v_s_tmp.num], s, strlen(s) + 1);
                v_s_tmp.num++;
                break;
            }
        }
        tmp++;
    }
    
    
    for (i = 0; i < v_s_tmp.num/2; i++) {
        snprintf(v_s_signed.val[i], sizeof(expression_t), "%s%s", v_s_tmp.val[2*i], v_s_tmp.val[2*i + 1]);
        v_s_signed.num++;
    }
    
    out->num = 0;
    
    for (i = 0; i < MAX_LEN; i++) {
        out->val[i].re = 0;
    }
    
    for (i = 0; i < v_s_signed.num; i++) {
        s = v_s_signed.val[i];
        if ((tmp = strstr(v_s_signed.val[i], "x^")) != NULL) {
            *tmp = 0;
            tmp += strlen("x^");
            degree = atof(tmp);
            if (degree > highest) {
                highest = degree;
            }
            if (*s == '-') {
                s++; out->val[degree].re = (*s == 0)?-1:-atof(s);
            } else if (*s == '+') {
                s++; out->val[degree].re = (*s == 0)?1:atof(s);
            }
        } else if ((tmp = strstr(v_s_signed.val[i], "x")) != NULL) {
            *tmp = 0;
            tmp += strlen("x");
            if (*s == '-') {
                s++; out->val[1].re = (*s == 0)?-1:-atof(s); highest = 1;
            } else if (*s == '+') {
                s++; out->val[1].re = (*s == 0)?1:atof(s); highest = 1;
            }
        } else {
            if (*s == '-') {
                s++; out->val[0].re = -atof(s);
            } else if (*s == '+') {
                s++; out->val[0].re = atof(s);
            }
        }
    }
    
    out->num = highest + 1;
    
    return 0;
}

int add_s(large_expression_t out, large_expression_t e1, large_expression_t e2)
{
    vector_t v1 = {0}, v2 = {0}, v3 = {0};
    unsigned int i, larger_num;
    
    arguments_vector(&v1, e1);
    arguments_vector(&v2, e2);
    
    if (v1.num >= v2.num) {
        larger_num = v1.num;
    } else {
        larger_num = v2.num;
    }
    
    v3.num = larger_num;
    
    for (i = 0; i < larger_num; i++) {
        v3.val[i].re = v1.val[i].re + v2.val[i].re;
    }
    
    arguments_to_polynomial_s(out, &v3);
    
    return 0;
}

int subtract_s(large_expression_t out, large_expression_t e1, large_expression_t e2)
{
    vector_t v1 = {0}, v2 = {0}, v3 = {0};
    unsigned int i, larger_num;
    
    arguments_vector(&v1, e1);
    arguments_vector(&v2, e2);
    
    if (v1.num >= v2.num) {
        larger_num = v1.num;
    } else {
        larger_num = v2.num;
    }
    
    v3.num = larger_num;
    
    for (i = 0; i < larger_num; i++) {
        v3.val[i].re = v1.val[i].re - v2.val[i].re;
    }
    
    arguments_to_polynomial_s(out, &v3);
    
    return 0;
}

int multiply_s(large_expression_t out, large_expression_t e1, large_expression_t e2)
{
    vector_t v1 = {0}, v2 = {0}, v3 = {0};
    unsigned int i;
    

    memset(out, 0, sizeof(large_expression_t));
    
    arguments_vector(&v1, e1);
    arguments_vector(&v2, e2);
    
    if ((v1.num > 1) && (v2.num > 1)) {
        multiply_polynomials(&v3, &v1, &v2);
    } else if ((v1.num == 1) && (v2.num > 1)) {
        for (i = 0; i < v2.num; i++) {
            v3.val[i].re = v1.val[0].re * v2.val[i].re;
            v3.num++;
        }
    } else if ((v1.num > 1) && (v2.num == 1)) {
        for (i = 0; i < v1.num; i++) {
            v3.val[i].re = v2.val[0].re * v1.val[i].re;
            v3.num++;
        }
    } else {
        v3.val[0].re = v2.val[0].re * v1.val[0].re;
        v3.num = 1;
    }
    
    arguments_to_polynomial_s(out, &v3);
    
   
    return 0;
}

int arguments_to_polynomial_s(large_expression_t out, vector_t *v)
{
    unsigned int i;
    expression_t tmp = {0};
    
    memset(out, 0, sizeof(large_expression_t));
    
    for (i = 0; i < v->num; i++) {
        if ((i == 1) && (fabs(v->val[i].re) != 0)) {
            if (v->val[i].re < 0) {
                snprintf(tmp, sizeof(expression_t), "- %0.4fx ", fabs(v->val[i].re));
            } else {
                snprintf(tmp, sizeof(expression_t), "+ %0.4fx ", v->val[i].re);
            }
            strncat(out, tmp, sizeof(expression_t));
        } else if ((i == 0) && (fabs(v->val[i].re) != 0)) {
            if (v->val[i].re < 0) {
                snprintf(tmp, sizeof(expression_t), "- %0.4f ", fabs(v->val[i].re));
            } else {
                snprintf(tmp, sizeof(expression_t), "+ %0.4f ", v->val[i].re);
            }
            strncat(out, tmp, sizeof(expression_t));
        } else if (fabs(v->val[i].re) != 0) {
            if (v->val[i].re < 0) {
                snprintf(tmp, sizeof(expression_t), "- %0.4fx^%d ", fabs(v->val[i].re), i);
            } else {
                snprintf(tmp, sizeof(expression_t), "+ %0.4fx^%d ", v->val[i].re, i);
            }
            strncat(out, tmp, sizeof(expression_t));
        }
        
    }
    
    return 0;
}

int determinant_s(large_expression_t out, matrix_s_t *in)
{
    unsigned int j;
    matrix_s_t m;
    large_expression_t temp, str, pr, pr1, pr2;

    if (in->rows != in->cols) {
        return -1;
    }
   
    if (in->rows == 1) {
        strncpy(out, in->val[0][0], strlen(in->val[0][0]) + 1);
        return 0;
    }

    if (in->rows == 2) {
        multiply_s(pr1, in->val[0][0], in->val[1][1]);
        multiply_s(pr2, in->val[0][1], in->val[1][0]);
        
        subtract_s(out, pr1, pr2);

        return 0;
    }
   
    *out = 0;

    for (j = 0; j < in->cols; j++) {
        minor_s(&m, in, 0, j);
        //print_matrix_s(&m);
        if (determinant_s(temp, &m) != 0) {
            wifi_util_error_print(WIFI_CSI,"%s:%d: Can not find determinant\n", __func__, __LINE__);
            return -1;
        }

        multiply_s(pr, in->val[0][j], temp);
        
        // out += out + pow(-1, j) * pr
        strncpy(str, out, strlen(out) + 1);
        
        if (pow(-1, j) > 0) {
            add_s(out, str, pr);
        } else {
            subtract_s(out, str, pr);
        }

        *temp = 0;
    }

    return 0;
}

int determinant(double *out, matrix_t *in)
{
    unsigned int j;
    double temp = 0.0;
    matrix_t m;

    if (in->rows != in->cols) {
        return -1;
    }

    if (in->rows == 1) {
        *out = in->val[0][0].re;
        return 0;
    }

    if (in->rows == 2) {
        *out = in->val[0][0].re * in->val[1][1].re - in->val[0][1].re * in->val[1][0].re;
        return 0;
    }

    *out = 0.0;

    for (j = 0; j < in->cols; j++) {
        minor(&m, in, 0, j);
        if (determinant(&temp, &m) != 0) {
            wifi_util_error_print(WIFI_CSI,"%s:%d: Can not find determinant\n", __func__, __LINE__);
            return -1;
        }
        *out += pow(-1, j) * in->val[0][j].re * temp;
        temp = 0.0;
    }

    return 0;
}

int cofactor(matrix_t *out, matrix_t *in)
{
    matrix_t m;
    unsigned int i, j, rows, cols;

    rows = in->rows;
    cols = in->cols;

    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            minor(&m, in, i, j);
            print_matrix(&m);
            //out->val[i][j] = pow(-1, i+j) * determinant(&m);
        }
    }

    out->rows = in->rows;
    out->cols = in->cols;

    return 0;
}

int polynomial(vector_t *out, vector_t *in)
{
    vector_t sum = {0}, temp = {0}, possible_root = {2};
    double root, lower_limit = -1000, d = 0;
    double upper_limit = 1000, increments = 0.0001;
    unsigned int i;

    if (in->num == 2) {
        out->val[out->num].re = (0.0 - in->val[1].re)/in->val[0].re;
        out->num++;
        return 0;
    }

    if (in->num == 3) {
        d = pow(in->val[1].re, 2) - 4*in->val[0].re*in->val[2].re;
        if (d < 0) {
            d *= -1;
            out->val[out->num].re = (-in->val[1].re)/(2*in->val[0].re);
            out->val[out->num].im = (sqrt(d))/(2*in->val[0].re);
            out->num++;
            out->val[out->num].re = (-in->val[1].re)/(2*in->val[0].re);
            out->val[out->num].im = -(sqrt(d))/(2*in->val[0].re);
            out->num++;
        } else {
            out->val[out->num].re = (-in->val[1].re + sqrt(d))/(2*in->val[0].re);
            out->num++;
            out->val[out->num].re = (-in->val[1].re - sqrt(d))/(2*in->val[0].re);
            out->num++;
        }
        
        return 0;
    }

    sum.num = in->num;
    sum.val[0].re = in->val[0].re;
    possible_root.val[1].re = 999;

    root = lower_limit;

    while (root < upper_limit) {

        for (i = 0; i < in->num - 1; i++) {
            sum.val[i + 1].re = in->val[i + 1].re + root*sum.val[i].re;
        }
        
        if (is_zero(sum.val[i].re, 1)) {
            break;
        } else if (is_zero(sum.val[i].re, 0)) {
            if (fabs(sum.val[i].re) < fabs(possible_root.val[1].re)) {
                possible_root.val[0].re = root;
                possible_root.val[1].re = sum.val[i].re;
            }
        }
        
        if (decimate(&root, 2) == 27.27) {
            wifi_util_error_print(WIFI_CSI,"root:%0.4f\tsum:%0.4f\n", root, sum.val[i].re);
        }

        root += increments;
    }

    if (root > upper_limit) {
        root = possible_root.val[0].re;
    }

    out->val[out->num].re = root;
    out->num++;

    temp.num = sum.num - 1;
    for (i = 0; i < temp.num; i++) {
        temp.val[i].re = root * sum.val[i].re;
    }

    temp.val[i - 1].re += sum.val[i].re;
    d = temp.val[0].re;
    for (i = 0; i < temp.num; i++) {
        temp.val[i].re /= d;
    }

    if (polynomial(out, &temp) != 0) {
        wifi_util_error_print(WIFI_CSI,"%s:%d:Cannot find roots\n", __func__, __LINE__);
        return -1;
    }

    return 0;
}

int quadratic(vector_t *out, vector_t *in)
{
    
    if (in->num != 3) {
        return -1;
    }

    out->num = 2;

    out->val[0].re = (-in->val[1].re + sqrt(pow(in->val[1].re, 2) - 4*in->val[0].re*in->val[2].re))/(2*in->val[0].re);
    out->val[1].re = (-in->val[1].re - sqrt(pow(in->val[1].re, 2) - 4*in->val[0].re*in->val[2].re))/(2*in->val[0].re);

    return 0;
}

int kurtosis(matrix_t *out, matrix_t *in)
{
    unsigned int i, j;
    matrix_t m1 = {0, 0}, mu4 = {0, 0}, mu2 = {0, 0};
    double mn;

    m1.rows = in->rows;
    m1.cols = 1;

    mu2.rows = 1;
    mu2.cols = in->cols;
    mu4.rows = 1;
    mu4.cols = in->cols;
    out->rows = 1;
    out->cols = in->cols;

    for (j = 0; j < in->cols; j++) {
        for (i = 0; i < in->rows; i++) {
            m1.val[i][0].re = in->val[i][j].re;
        }
        if (mean(&mn, &m1) != 0) {
            wifi_util_error_print(WIFI_CSI,"%s:%d: Cannot calculate mean\n", __func__, __LINE__);
            return -1;
        }

        for (i = 0; i < in->rows; i++) {
            mu4.val[0][j].re += pow((in->val[i][j].re - mn), 4);
            mu2.val[0][j].re += pow((in->val[i][j].re - mn), 2);
        }

        mu4.val[0][j].re = mu4.val[0][j].re/in->rows;
        mu2.val[0][j].re = mu2.val[0][j].re/in->rows;

        out->val[0][j].re = mu4.val[0][j].re/pow(mu2.val[0][j].re, 2);

    }

    return 0;
}

int correlation(matrix_t *out, matrix_t *in)
{
    unsigned int i, j;
    matrix_t m1 = {0, 0}, m2 = {0, 0}, m3 = {0, 0}, m4 = {0, 0};
    double m, s;

    m2.rows = in->rows;
    m2.cols = in->cols;

    m1.rows = in->rows;
    m1.cols = 1;

    for (j = 0; j < in->cols; j++) {
        for (i = 0; i < in->rows; i++) {
            m1.val[i][0] = in->val[i][j];
        }
        if (mean(&m, &m1) != 0) {
            wifi_util_error_print(WIFI_CSI,"%s:%d: Cannot calculate mean\n", __func__, __LINE__);
            return -1;
        }

        if (stddev(&s, &m1) != 0) {
            wifi_util_error_print(WIFI_CSI,"%s:%d: Cannot calculate standard deviation\n", __func__, __LINE__);
            return -1;
        }

        for (i = 0; i < in->rows; i++) {
            m2.val[i][j].re = (in->val[i][j].re - m)/s;
        }

    }

    if (transpose(&m3, &m2) != 0) {
        wifi_util_error_print(WIFI_CSI,"%s:%d: Cannot calculate transpose\n", __func__, __LINE__);
        return -1;
    }
   
    if (multiply(&m4, &m3, &m2) != 0) {
        wifi_util_error_print(WIFI_CSI,"%s:%d: Cannot multiply matrices\n", __func__, __LINE__);
        return -1;
    }

    out->rows = m4.rows;
    out->cols = m4.cols;

    for (i = 0; i < in->rows; i++) {
        for (j = 0; j < in->cols; j++) {
            out->val[i][j].re = m4.val[i][j].re/(in->rows - 1);
        }
    }

    return 0;
}

int covariance(matrix_t *out, matrix_t *in)
{
    unsigned int i, j;
    double m;
    matrix_t m1 = {0, 0}, m2 = {0, 0}, m3 = {0, 0}, m4 = {0, 0};

    if (in->rows == 1) {
        wifi_util_error_print(WIFI_CSI,"%s:%d: Cannot calculate covariance of insufficient valued matrices\n", __func__, __LINE__);
        return -1;
    }

    m2.rows = in->rows;
    m2.cols = in->cols;

    m1.rows = m2.rows;
    m1.cols = 1;

    for (j = 0; j < in->cols; j++) {
        for (i = 0; i < m1.rows; i++) {
            m1.val[i][j] = in->val[i][j];
        }
   
        if (mean(&m, &m1) != 0) {
            wifi_util_error_print(WIFI_CSI,"%s:%d: Cannot calculate mean\n", __func__, __LINE__);
            return -1;
        }

        for (i = 0; i < in->rows; i++) {
            m2.val[i][j].re = in->val[i][j].re - m;
        }

    }
   
    if (transpose(&m3, &m2) != 0) {
        wifi_util_error_print(WIFI_CSI,"%s:%d: Cannot calculate transpose\n", __func__, __LINE__);
        return -1;
    }
   
    if (multiply(&m4, &m3, &m2) != 0) {
        wifi_util_error_print(WIFI_CSI,"%s:%d: Cannot multiply matrices\n", __func__, __LINE__);
        return -1;
    }


    out->rows = m4.rows;
    out->cols = m4.cols;

    for (i = 0; i < m3.rows; i++) {
        for (j = 0; j < m3.cols; j++) {
            out->val[i][j].re = m3.val[i][j].re/(in->rows - 1);
        }
    }

    return 0;
}

int stddev(double *out, matrix_t *in)
{
    double d;

    if (variance(&d, in) != 0) {
        wifi_util_error_print(WIFI_CSI,"%s:%d: Cannot calculate variance\n", __func__, __LINE__);
        return -1;
    }

    *out = sqrt(d);

    return 0;
}

int variance(double *out, matrix_t *in)
{
    unsigned int i;
    matrix_t m1 = {0, 0}, m2 = {0, 0}, m3 = {0, 0}, m4 = {0, 0};
    double m;

    if (in->rows == 1) {
        *out = in->val[0][0].re;
        return 0;
    }

    m1.rows = in->rows;
    m1.cols = 1;

    for (i = 0; i < m1.rows; i++) {
        m1.val[i][0].re = in->val[i][0].re;
    }

    if (mean(&m, &m1) != 0) {
        wifi_util_error_print(WIFI_CSI,"%s:%d: Cannot calculate mean\n", __func__, __LINE__);
        return -1;
    }

    m2.cols = 1;
    m2.rows = in->rows;

    for (i = 0; i < in->rows; i++) {
        m2.val[i][0].re = in->val[i][0].re - m;
    }

    if (transpose(&m3, &m2) != 0) {
        wifi_util_error_print(WIFI_CSI,"%s:%d: Cannot calculate transpose\n", __func__, __LINE__);
        return -1;
    }
   
    if (multiply(&m4, &m3, &m2) != 0) {
        wifi_util_error_print(WIFI_CSI,"%s:%d: Cannot multiply matrices\n", __func__, __LINE__);
        return -1;
    }

    *out = m4.val[0][0].re/(in->rows - 1);

    return 0;
}

int mean(double *out, matrix_t *in)
{
    double sum = 0;
    unsigned int i;

    if ((in->rows == 0) || (in->cols == 0)) {
        return 0;
    }

    for (i = 0; i < in->rows; i++) {
        sum += in->val[i][0].re;
    }

    *out = sum/in->rows;
    return 0;
}

int transpose(matrix_t *out, matrix_t *in)
{
    unsigned int i = 0, j = 0;

    out->rows = in->cols;
    out->cols = in->rows;

    for (i = 0; i < out->rows; i++) {
        for (j = 0; j < out->cols; j++) {
            out->val[i][j] = in->val[j][i];
        }
    }
    return 0;
}

int multiply(matrix_t *out, matrix_t *m1, matrix_t *m2)
{
    unsigned int i = 0, j = 0, k = 0;

    if (m1->cols != m2->rows) {
        wifi_util_error_print(WIFI_CSI,"%s:%d: Matrices can't be multipled, mismatch of 1st Matrix Columns: %d and 2nd Matrix Rows: %d\n", __func__, __LINE__,
                    m1->cols, m2->rows);
        return -1;
    }

    out->rows = m1->rows;
    out->cols = m2->cols;

    for (i = 0; i < out->rows; i++) {
        for (j = 0; j < out->cols; j++) {
            for (k = 0; k < m1->cols; k++) {
                out->val[i][j].re += m1->val[i][k].re * m2->val[k][j].re;
            }
        }
    }

    return 0;
}

int linear_eq(matrix_t *out, matrix_t *a, matrix_t *b)
{
    matrix_t inv;
    
    if (inverse(&inv, a) != 0) {
        wifi_util_error_print(WIFI_CSI,"%s:%d: Linear equation can not be solved\n", __func__, __LINE__);
        return -1;
    }
    
    multiply(out, &inv, b);
    
    return 0;
}

int eigen_values(vector_t *eigen, matrix_t *in)
{
    matrix_s_t m;
    large_expression_t e;
    unsigned int i, j;
    double d;
    vector_t args = {0};
    
    m.rows = in->rows;
    m.cols = in->cols;
    for (i = 0; i < in->rows; i++) {
        for (j = 0; j < in->cols; j++) {
            if (i == j) {
                snprintf(m.val[i][j], sizeof(expression_t), "%0.4f - x", in->val[i][j].re);
            } else {
                snprintf(m.val[i][j], sizeof(expression_t), "%0.4f", in->val[i][j].re);
            }
        }
    }
    
    if (determinant_s(e, &m) != 0) {
        wifi_util_error_print(WIFI_CSI,"%s:%d: Can not find determinant\n", __func__, __LINE__);
        return -1;
    }
    
    arguments_vector(&args, e);
    for (i = 0; i < args.num/2 + 1; i++) {
        d = args.val[args.num - i - 1].re;
        args.val[args.num - i - 1].re = args.val[i].re;
        args.val[i].re = d;
    }
    
    print_vector(&args);
    
    if (polynomial(eigen, &args) != 0) {
        wifi_util_error_print(WIFI_CSI,"%s:%d: Can not find eigen values\n", __func__, __LINE__);
        return -1;
    }
    
    return 0;
}

int eigen_vectors(vector_t *eigen_vec, number_t *eigen_val, matrix_t *in)
{
    matrix_t m, n, null, out;
    unsigned int i, j;
    
    m.rows = in->rows;
    m.cols = in->cols;
    for (i = 0; i < in->rows; i++) {
        for (j = 0; j < in->cols; j++) {
            if (i == j) {
                m.val[i][j].re = in->val[i][j].re - eigen_val->re;
            } else {
                m.val[i][j].re = in->val[i][j].re;
            }
        }
    }
    
    print_matrix(&m);
    if (inverse(&n, &m) != 0) {
        wifi_util_error_print(WIFI_CSI,"%s:%d: Cannot inverse matrix\n", __func__, __LINE__);
        return -1;
    }
    wifi_util_info_print(WIFI_CSI,"Inversed matrix\n");
    print_matrix(&n);
    
    null.cols = 1;
    null.rows = n.rows;
    
    multiply(&out, &n, &null);
    
    wifi_util_info_print(WIFI_CSI,"Output\n");
    print_matrix(&out);
    
    return 0;
}

int eigen_decompose(vector_t *eigen, vector_array_t *eigen_vecs, matrix_t *in)
{
    unsigned int i;
    vector_t eigen_vec;
    
    for (i = 0; i < MAX_LEN; i++) {
        memset(&eigen->val[i], 0, sizeof(number_t));
    }

    if (in == NULL) {
        return -1;
    }

    if (in->rows != in->cols) {
        return -1;
    }

    if (in->rows == 0 || in->rows >= MAX_LEN) {
        return -1;
    }
    
    if (eigen_values(eigen, in) != 0) {
        wifi_util_error_print(WIFI_CSI,"%s:%d: Can not find eigen values\n", __func__, __LINE__);
        return -1;
    }
    
    //print_vector(eigen);
    
    for (i = 0; i < eigen->num; i++) {
        if ((eigen->val[i].im == 0) && (eigen_vectors(&eigen_vec, &eigen->val[i], in) != 0)) {
            
        }
    }
    
    return 0;
}

int s_to_vector(vector_t *out, vector_s_t *in)
{
    unsigned int i;
    
    
    for (i = 0; i < in->num; i++) {
        
        out->val[i].re = atof(in->val[i]);
        out->num++;
       
    }
    
    return 0;
}

int vector_to_s(vector_s_t *out, vector_t *in)
{
    unsigned int i = 0;
    
    for (i = 0; i < in->num; i++) {
        snprintf(out->val[out->num++], sizeof(expression_t), "%0.4f", in->val[i].re);
    }
    
    return 0;
}

int vector_s(vector_s_t *v, ...)
{
    va_list list;
    unsigned int i;
    char *s;

    if (v == NULL) {
        return -1;
    }

    if ((v->num == 0) || (v->num >= MAX_LEN)) {
        return -1;
    }

    va_start(list, v);

    for (i = 0; i < v->num; i++) {
        s = va_arg(list, char *);
        strncpy(v->val[i], s, strlen(s) + 1);
    }

    va_end(list);

    return 0;
}

int vector(vector_t *v, ...)
{
    va_list list;
    unsigned int i;

    if (v == NULL) {
        return -1;
    }

    if ((v->num == 0) || (v->num >= MAX_LEN)) {
        return -1;
    }

    va_start(list, v);

    for (i = 0; i < v->num; i++) {
        v->val[i].re  = va_arg(list, double);
    }

    va_end(list);

    return 0;
}

int matrix_s(matrix_s_t *m, ...)
{
    va_list list;
    unsigned int i, j;
    char *s;

    if (m == NULL) {
        return -1;
    }

    if ((m->rows == 0) || (m->rows >= MAX_LEN)) {
        return -1;
    }

    if ((m->cols == 0) || (m->cols >= MAX_LEN)) {
        return -1;
    }

    va_start(list, m);

    for (i = 0; i < m->rows; i++) {
        for (j = 0; j < m->cols; j++) {
            s = va_arg(list, char *);
            strncpy(m->val[i][j], s, strlen(s) + 1);
        }
    }
   
    va_end(list);

    return 0;
}

int matrix(matrix_t *m, ...)
{
    va_list list;
    unsigned int i, j;

    if (m == NULL) {
        return -1;
    }

    if ((m->rows == 0) || (m->rows >= MAX_LEN)) {
        return -1;
    }

    if ((m->cols == 0) || (m->cols >= MAX_LEN)) {
        return -1;
    }

    va_start(list, m);

    for (i = 0; i < m->rows; i++) {
        for (j = 0; j < m->cols; j++) {
            m->val[i][j].re = va_arg(list, double);
        }
    }
   
    va_end(list);

    return 0;
}

void print_number(number_t *n)
{
    if ((is_zero(n->re, 4) == true) && (is_zero(n->im, 4) == false)) {
        wifi_util_info_print(WIFI_CSI,"%0.4fi\t", n->im);
    } else if ((is_zero(n->im, 4) == true) && (is_zero(n->re, 4) == false)) {
        wifi_util_info_print(WIFI_CSI,"%0.4f\t", n->re);
    } else if ((is_zero(n->im, 4) == false) && (is_zero(n->re, 4) == false)) {
        if (n->im < 0) {
            wifi_util_info_print(WIFI_CSI,"%0.4f%0.4fi\t", n->re, n->im);
        } else {
            wifi_util_info_print(WIFI_CSI,"%0.4f+%0.4fi\t", n->re, n->im);
        }
    } else {
        wifi_util_info_print(WIFI_CSI,"0.0000\t");
    }
}

void print_vector(vector_t *v)
{
    unsigned int i;

    wifi_util_info_print(WIFI_CSI,"Vector:");
    for (i = 0; i < v->num; i++) {
        print_number(&v->val[i]);
    }

    wifi_util_info_print(WIFI_CSI,"\n");
}

void print_vector_s(vector_s_t *v)
{
    unsigned int i;

    wifi_util_info_print(WIFI_CSI,"Vector:");
    for (i = 0; i < v->num; i++) {
        wifi_util_info_print(WIFI_CSI,"%s\t", v->val[i]);
    }

    wifi_util_info_print(WIFI_CSI,"\n");
}

void print_matrix_s(matrix_s_t *m)
{
    unsigned int i, j;

    wifi_util_info_print(WIFI_CSI,"Matrix:");
    wifi_util_info_print(WIFI_CSI,"\n");
    for (i = 0; i < m->rows; i++) {
        for (j = 0; j < m->cols; j++) {
            wifi_util_info_print(WIFI_CSI,"%s\t", m->val[i][j]);
        }
        wifi_util_info_print(WIFI_CSI,"\n");
    }

    wifi_util_info_print(WIFI_CSI,"\n");

}

void print_matrix(matrix_t *m)
{
    unsigned int i, j;

    wifi_util_info_print(WIFI_CSI,"Matrix:\r\n");
    for (i = 0; i < m->rows; i++) {
        for (j = 0; j < m->cols; j++) {
            print_number(&m->val[i][j]);
        }
        wifi_util_info_print(WIFI_CSI,"\n");
    }

    wifi_util_info_print(WIFI_CSI,"\n");
}

bool is_zero(double x, int n) {
    if (n == 0) {
        return abs((int)x) == 0;
    }
    
    return fabs(x) < pow(10, -n);
}

double decimate(double *in, int pos)
{
    double d = (*in) * pow(10, pos);
    
    return abs((int)d) * pow(10, -pos);
}
