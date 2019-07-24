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


 */

#include "basis.h"
#include <math.h>

double SIntegral(double a[3],double b[3],int,int,int,int,int,int,double,double);
double JIntegral(double a[3],double b[3],int,int,int,int,int,int,double,double);

double orbital_SIntegral(orbital,orbital);
double orbital_JIntegral(orbital,orbital);
double orbital_GIntegral(orbital,orbital);
double orbital_FIntegral(orbital,orbital);
// double orbital_GJIntegral(orbital,orbital);
// double orbital_FJIntegral(orbital,orbital);
// double orbital_FFIntegral(orbital,orbital);
