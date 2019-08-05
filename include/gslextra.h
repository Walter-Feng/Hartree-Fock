#ifndef __GSLEXTRA_H__
#define __GSLEXTRA_H__

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>

#endif


void gsl_vector_complex_convert(gsl_vector * source, gsl_vector_complex * target, int length);
void gsl_matrix_complex_convert(gsl_matrix * source, gsl_matrix_complex * target, int rows, int columns);
void gsl_vector_complex_extract(gsl_vector_complex * source, gsl_vector * real, gsl_vector * imag, int length);
void gsl_matrix_complex_extract(gsl_vector_complex * source, gsl_matrix * real,gsl_matrix * imag, int rows, int columns);
void gsl_vector_complex_combine(gsl_vector * real, gsl_vector * imag, gsl_vector_complex * target);
void gsl_matrix_complex_combine(gsl_matrix * real,gsl_matrix * imag, gsl_matrix_complex * target);
void gsl_matrix_diag(gsl_matrix * target, gsl_vector * diag, int length);
void gsl_matrix_complex_diag(gsl_matrix_complex * target, gsl_vector_complex * diag, int length);
void gsl_matrix_mul(gsl_matrix * A, gsl_matrix *B, gsl_matrix * Result,int Acolumn,int Arow,int Bcolumn);
void gsl_matrix_complex_mul(gsl_matrix_complex * A, gsl_matrix_complex *B, gsl_matrix_complex * Result,int Acolumn,int Arow,int Bcolumn);
double gsl_vector_inner_product(gsl_vector * A, gsl_vector * B,int length);
gsl_complex gsl_vector_complex_product(gsl_vector_complex * A, gsl_vector_complex * B, int length);
gsl_complex gsl_vector_complex_inner_product(gsl_vector_complex * A, gsl_vector_complex * B, int length);
void gsl_vector_transform(gsl_vector * vec,gsl_matrix * trf,int length);
void gsl_vector_complex_transform(gsl_vector_complex * vec, gsl_matrix_complex * trf,int length);
void gsl_matrix_unitmatrix(gsl_matrix * m,int length);
void gsl_matrix_complex_unitmatrix(gsl_matrix_complex * m,int length);
void gsl_vector_complex_conjugate(gsl_vector_complex * v, int length);
void gsl_matrix_complex_conjugate(gsl_matrix_complex * m, int rows, int columns);

void gsl_eigen_Lowdin_diag(gsl_matrix * m, gsl_matrix * S, gsl_vector * eigen, gsl_matrix * eigenvec);