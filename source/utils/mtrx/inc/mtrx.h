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

/*
 * NOTE: This file is included also in OneWifi project which is C based, thus
 * there should be only usage of C based constructs in this file.
 * C++ constructs are not allowed in this file.
 */
#ifndef MTRX_H
#define MTRX_H

#ifdef __cplusplus
extern "C"
{
#endif

#define     MAX_LEN     64

typedef char expression_t[32];
typedef char large_expression_t[2048];

typedef struct {
    double re;
    double im;
} number_t;

typedef struct {
    unsigned int num;
    number_t  val[MAX_LEN];
} vector_t;

typedef struct {
    unsigned int num;
    expression_t    val[MAX_LEN];
} vector_s_t;

typedef struct {
    unsigned int    rows;
    unsigned int    cols;
    number_t  val[MAX_LEN][MAX_LEN];
} matrix_t;

typedef struct {
    unsigned int    rows;
    unsigned int    cols;
    expression_t    val[MAX_LEN][MAX_LEN];
} matrix_s_t;

typedef matrix_t vector_array_t;

bool is_zero(double x, int n);
int vector(vector_t *v, ...);
int vector_s(vector_s_t *v, ...);
int vector_to_s(vector_s_t *out, vector_t *in);
int s_to_vector(vector_t *out, vector_s_t *in);
void print_number(number_t *n);
void print_vector(vector_t *v);
void print_vector_s(vector_s_t *v);
void print_matrix(matrix_t *m);
void print_matrix_s(matrix_s_t *m);
int matrix(matrix_t *m, ...);
int matrix_s(matrix_s_t *m, ...);
int arguments_vector(vector_t *v, large_expression_t e);
int multiply(matrix_t *out, matrix_t *m1, matrix_t *m2);
int multiply_s(large_expression_t out, large_expression_t e1, large_expression_t e2);
int subtract_s(large_expression_t out, large_expression_t e1, large_expression_t e2);
int add_s(large_expression_t out, large_expression_t e1, large_expression_t e2);
int arguments_to_polynomial_s(large_expression_t out, vector_t *v);
int transpose(matrix_t *out, matrix_t *in);
int mean(double *out, matrix_t *in);
int variance(double *out, matrix_t *in);
int stddev(double *out, matrix_t *in);
int covariance(matrix_t *out, matrix_t *in);
int correlation(matrix_t *out, matrix_t *in);
int kurtosis(matrix_t *out, matrix_t *in);
int quadratic(vector_t *out, vector_t *in);
int polynomial(vector_t *out, vector_t *in);
int cofactor(matrix_t *in, matrix_t *out);
int minor(matrix_t *out, matrix_t *in, unsigned int row, unsigned int col);
int minor_s(matrix_s_t *out, matrix_s_t *in, unsigned int row, unsigned int col);
int determinant(double *out, matrix_t *in);
int determinant_s(large_expression_t out, matrix_s_t *in);
int clone(matrix_t *out, matrix_t *in);
int adjoint(matrix_t *out, matrix_t *in);
int inverse(matrix_t *out, matrix_t *in);
int multiply_polynomials(vector_t *out, vector_t *v1, vector_t *v2);
int linear_eq(matrix_t *out, matrix_t *a, matrix_t *b);
double decimate(double *in, int pos);

int nullify(matrix_t *in);
int nullify_s(matrix_s_t *in);

int eigen_values(vector_t *eigen, matrix_t *in);
int eigen_vector(vector_t *eigen_vec, number_t *eigen_val, matrix_t *in);
int eigen_decompose(vector_t *eigen, vector_array_t *eigen_vecs, matrix_t *in);

#ifdef __cplusplus
}
#endif

#endif // MTRX_H
