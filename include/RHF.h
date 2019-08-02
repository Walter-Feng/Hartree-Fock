#ifndef __RHF_H__
#define __RHF_H__
#endif

#include "integral.h"

#include <gsl/gsl_eigen.h>


double fock_matrix_element(orbital * a, orbital * b, orbital * HEAD, gsl_matrix * coef, int length, int el_num);

void fock_matrix(gsl_matrix * dest, gsl_matrix * coef, orbital * HEAD, int length, int el_num);
