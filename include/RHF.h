#ifndef __RHF_H__
#define __RHF_H__

#include "integral.h"
#include "../include/basis.h"
#include "../include/gslextra.h"
#include "../include/gslprint.h"
#include <gsl/gsl_eigen.h>

#endif


double fock_matrix_element(orbital * a, orbital * b, orbital * HEAD, gsl_matrix * coef, int length, int el_num);

void fock_matrix(gsl_matrix * dest, gsl_matrix * coef, orbital * HEAD, int length, int el_num);

int RHF_SCF_print(gsl_vector * energy, gsl_matrix * coef, orbital * HEAD, int length, int el_num, int iteration_max, double errmax, int countmax);
