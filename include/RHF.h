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

double nuclear_attraction_energy_matrix_element(orbital * a, orbital * b, atomic_orbital * atom_HEAD);

double single_electron_hamiltonian_matrix_element(orbital * a, orbital * b, atomic_orbital * atom_HEAD);

double fock_matrix_element(orbital * a, orbital * b, orbital * HEAD, atomic_orbital * atom_HEAD, gsl_matrix * coef, int length, int el_num);

void kinetic_energy_matrix(gsl_matrix * dest, orbital * HEAD, int length);

void nuclear_attraction_energy_matrix(gsl_matrix * dest, orbital * HEAD, atomic_orbital * atom_HEAD, int length);

void core_hamiltonian_matrix(gsl_matrix * dest, orbital * HEAD, atomic_orbital * atom_HEAD, int length);

void two_electron_quad_tensor(gsl_quad_tensor * dest, orbital * HEAD, int length);

void fock_matrix(gsl_matrix * dest, gsl_matrix * coef, orbital * HEAD, atomic_orbital * atom_HEAD, int length, int el_num);

void initial_guess(gsl_matrix * dest, gsl_matrix * S, orbital * HEAD, atomic_orbital * atom_HEAD, int length);

int RHF_SCF_print(gsl_vector * energy, gsl_matrix * coef, orbital * HEAD, atomic_orbital * atom_HEAD, int length, int el_num, int iteration_max, double errmax, int countmax, double alpha);

double nuclei_repulsion(atomic_orbital * atomlist_HEAD);