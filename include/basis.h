#ifndef __BASIS_H__
#define __BASIS_H__
#endif

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <string.h>

//This struct is aimed to store each electronic shell of an atom
typedef struct orbital
{
    //Angular quantum number;
    int L;
    //magnetic quantum number;
    int m;
    //main quantum number;
    int n;
    //the name of the orbital;
    char label[10];

    //the angular momentum exponents combined with their coefficients (concerning the normalization for a linear combination of terms with different angular momentum exponents)
    struct  angcoef{
        int a[3];
        double coef;
    }A[4];

    //This will tell the functions how long is the angcoef
    int length;

    //the total number of the terms of coefficients and exponents 
    int total;

    //The list of exponents
    double * exponents;
    //THe corresponding coefficient list, which will be normalized during the process;
    double * coefficients;

    //The cartesian coordinate of the center of the orbital;
    double cartesian[3];

    //pointer storing the next orbital
    orbital* NEXT;
}orbital;

typedef struct atomic_orbital
{
    // The atomic number of the atom;
    int N;
    // The name of the atom;
    char name[5];

    //The cartesian coordinate of the atom;
    double cartesian[3];

    //The HEAD of the orbital
    orbital * orbital_HEAD;

    //Next atom
    atomic_orbital * NEXT;
}atomic_orbital;

orbital* orbital_calloc(size_t);
void orbital_free(orbital *);

atomic_orbital* atomic_orbital_calloc();
void atomic_orbital_free(atomic_orbital *);

void basis_fscanf(FILE *,atomic_orbital *);

void orbital_label(char *,int,int,int);
void orbital_angcoef_set(orbital *);

void atomic_orbital_cpy(atomic_orbital *, atomic_orbital *);
void orbital_cpy(orbital *, orbital *);

int orbital_count(orbital *);