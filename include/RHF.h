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

double fock_matrix_element(gsl_quad_tensor * v, gsl_matrix * density_matrix, gsl_matrix * h_matrix, int i, int j, int length);

double HF_energy(gsl_quad_tensor * v, gsl_matrix * density_matrix, gsl_matrix * h_matrix, int length);

void kinetic_energy_matrix(gsl_matrix * dest, orbital * HEAD, int length);

void nuclear_attraction_energy_matrix(gsl_matrix * dest, orbital * HEAD, atomic_orbital * atom_HEAD, int length);

void core_hamiltonian_matrix(gsl_matrix * dest, orbital * HEAD, atomic_orbital * atom_HEAD, int length);

//obtain the quad tensor of the two-electron Coulomb integrals 
void two_electron_quad_tensor(gsl_quad_tensor * dest, orbital * HEAD, int length);

//obtain the fock matrix
void fock_matrix(gsl_matrix * dest, gsl_quad_tensor * v, gsl_matrix * density_matrix, gsl_matrix * h_matrix, int length);

// perform initial guess of the coefficient matrix by performing diagonalization of core hamiltonian matrix
void initial_guess(gsl_matrix * dest, gsl_matrix * core_hamiltonian, gsl_matrix * S, int length);

// calculate density matrix from coefficient matrix and the number of electrons
void density_matrix(gsl_matrix * dest, gsl_matrix * coef, int el_num, int length);

//perform RHF as well as printing information
int RHF_SCF_print(double * tot_energy, gsl_vector * energy, gsl_matrix * coef, orbital * HEAD, atomic_orbital * atom_HEAD, int length, int el_num, int iteration_max, double errmax, int countmax, double alpha, int mixing_type, int SCF_INITIAL_FLAG, int SCF_FOCK_FLAG, int SCF_COEF_FLAG, int FOCK_FLAG);

//calculate the nuclei repulsion energy
double nuclei_repulsion(atomic_orbital * atomlist_HEAD);