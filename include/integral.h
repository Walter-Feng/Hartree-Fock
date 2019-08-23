/*
    This section is aimed to calculate the gaussian integrals that occurs in everywhere

    This part is greatly inspired by May, Andrew James. Density fitting in explicitly correlated electronic structure theory. Diss. University of Bristol, 2006.
 */

/*

Here is the table of the relation between the names and the definition:

S   <a|b>
J   (a|b) = <a|1/r_{12}|b>
G   (a|e^{-\gamma r_{12}^2}|b)
F   (a|f_{12}|b) = (a|c_i exp(-\gamma_i r_{12}^2)|b)
GJ  (a|e^{-\gamma r_{12}^2 r_12^{-1}}|b)
FJ  (a|f_{12} r_{12}^{-1}|b)
FF  (a|f_{12}^2|b)
FT  (ab|[\hat{t}_1, f_{12}]|c) = (ab|[- \frac{1}{2} \nabla ^2, f_{12}]|c)
FTF (a|\frac{1}{2}[f_12,[\hat{t}_1 + \hat{t}_2, f_{12}]]|b)
F-F (a|f_{12}|b|f_{23}|c)
J-F (a|r_{12}^{-1}|b|f_{23}|c)
X   (ab|[\hat{t}_1,r_{12}]|c)
Y   (ab|[\hat{t}_1,r_{12}^{-1}])
Z   (a|\frac{1}{r_{1Z}}|b)

Currently only S integrals and J integrals are supported.

 */

#ifndef __BASIS_H__
#define __BASIS_H__
#include "basis.h"
#endif

#include <math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf.h>

typedef struct gaussian_chain{
    double R[3];
    int a[3];
    
    double exponent;

    double coefficient;

    gaussian_chain * NEXT;

}gaussian_chain;

gaussian_chain * gaussian_chain_calloc();

void gaussian_chain_free(gaussian_chain * HEAD);

double gaussian_chain_get(gaussian_chain * HEAD,double x, double y, double z);
double orbital_get(orbital * orbital,double x, double y, double z);

double Gamma(double z);
double Boys(double x, int n);
double Binomials(int n, int k);

double f(int k, int a, int b, double PA, double PB);
double tranformation_coefficient(int a[3], int b[3], int p[3], double PA[3], double PB[3], double xi, double AB);

double SIntegral(double ra[3], double rb[3], int ax, int ay, int az, int bx, int by, int bz, double alpha,double beta);
double JIntegral(double ra[3], double rb[3], int ax, int ay, int az, int bx, int by, int bz, double alpha,double beta, int m);
double ZIntegral(double ra[3], double rb[3], double rz[3], int ax, int ay, int az, int bx, int by, int bz, double alpha, double beta, int m);

double gaussian_chain_SIntegral(gaussian_chain * a, gaussian_chain * b);
double gaussian_chain_JIntegral(gaussian_chain * a, gaussian_chain * b);
double gaussian_chain_ZIntegral(gaussian_chain * a, gaussian_chain * b, double rz[3]);
double gaussian_chain_full_SIntegral(gaussian_chain * a_HEAD, gaussian_chain * b_HEAD);
double gaussian_chain_full_JIntegral(gaussian_chain * a_HEAD, gaussian_chain * b_HEAD);
double gaussian_chain_full_ZIntegral(gaussian_chain * a_HEAD, gaussian_chain * b_HEAD, double rz[3]);

void gaussian_chain_derivative(gaussian_chain * dest, gaussian_chain * src, int key);
void gaussian_chain_second_derivative(gaussian_chain * dest, gaussian_chain * src, int key);
void gaussian_chain_laplacian(gaussian_chain * dest, gaussian_chain * src);

double gaussian_chain_kinetic_energy(gaussian_chain * a_HEAD, gaussian_chain * b_HEAD);

void single_electron_transform(gaussian_chain * HEAD, orbital * a);
void two_electron_transform(gaussian_chain * HEAD, orbital * a, orbital * b);

double orbital_SIntegral(orbital * a,orbital * b);
double orbital_JIntegral(orbital * a,orbital * b);
double orbital_ZIntegral(orbital * a,orbital * b, double rz[3]);
//double orbital_GIntegral(orbital *,orbital *);
//double orbital_FIntegral(orbital *,orbital *);
//double orbital_GJIntegral(orbital *,orbital *);
//double orbital_FJIntegral(orbital *,orbital *);
//double orbital_FFIntegral(orbital *,orbital *);

double two_electron_JIntegral(orbital * a_1, orbital * b_1, orbital * c_2, orbital * d_2);

double orbital_kinetic_energy(orbital * a, orbital * b);

void orbital_S_matrix(gsl_matrix * dest, orbital * HEAD);