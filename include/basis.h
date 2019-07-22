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
    //magnetic quantum numberl
    int m;
    //main quantum number;
    int n;
    //the name of the orbital;
    char* label;

    int total;

    //The list of exponents
    double * exponents;
    //THe corresponding coefficient list, which will be normalized during the process;
    double * coefficients;

    //pointer storing the next orbital
    void *NEXT;
}orbital;

typedef struct atomic_orbital
{
    // The atomic number of the atom;
    int N;
    // The name of the atom;
    char* name;

    //The cartesian coordinate of the atom;
    double cartesian[3];

    //The HEAD of the orbital
    void* orbital_HEAD;

    //Next atom
    void* NEXT;
}atomic_orbital;

orbital* orbital_calloc(size_t);
void orbital_free(orbital *);

atomic_orbital* atomic_orbital_calloc();
void atomic_orbital_free(atomic_orbital *);

void basis_fscanf(FILE *,atomic_orbital *);

