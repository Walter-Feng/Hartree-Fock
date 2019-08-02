#ifndef __GSLPRINT_H__
#define __GSLPRINT_H__

#include <stdio.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_integration.h>

#endif

void gsl_matrix_fprint(FILE * stream, gsl_matrix * target,int rows,int columns, char const * format);
void gsl_vector_fprint(FILE * stream, gsl_vector * target,int length, char const * format);
void gsl_matrix_printf(gsl_matrix * m, int rows, int columns, char const * format);
void gsl_vector_printf(gsl_vector * m,int length, char const * format);
void gsl_vector_complex_fprint(FILE * stream, gsl_vector_complex * target, int length, char const * format);
void gsl_matrix_complex_fprint(FILE * stream, gsl_matrix_complex * target, int rows, int columns, char const * format);
void gsl_matrix_complex_printf(gsl_matrix_complex * target, int rows, int columns, char const * format);
void gsl_vector_complex_fprint(gsl_vector_complex * target, int length, char const * format);
