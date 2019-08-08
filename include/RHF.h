#ifndef __INTEGRAL_H__
#define __INTEGRAL_H__
#include "integral.h"
#endif

#ifndef __GSLEXTRA_H__
#define __GSLEXTRA_H__
#include "../include/gslextra.h"
#endif

#ifndef __GSLPRINT_H__
#define __GSLPRINT_H__
#include "../include/gslprint.h"
#endif

#include <gsl/gsl_eigen.h>


double fock_matrix_element(orbital * a, orbital * b, orbital * HEAD, atomic_orbital * atom_HEAD, gsl_matrix * coef, int length, int el_num);

void fock_matrix(gsl_matrix * dest, gsl_matrix * coef, orbital * HEAD, atomic_orbital * atom_HEAD, int length, int el_num);

int RHF_SCF_print(gsl_vector * energy, gsl_matrix * coef, orbital * HEAD, atomic_orbital * atom_HEAD, int length, int el_num, int iteration_max, double errmax, int countmax);

double nuclei_repulsion(atomic_orbital * atomlist_HEAD);