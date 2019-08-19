#define __BASIS_H__

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
    char label[20];

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
    //The corresponding coefficient list, which will be normalized during the process;
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

//To allocate the memory for struct orbital
orbital* orbital_calloc(int);
// free the memory for struct orbital
void orbital_free(orbital * HEAD);

orbital * orbital_enquiry(orbital * HEAD, int index);

//To allocate the memory for struct atomic_orbital
atomic_orbital* atomic_orbital_calloc();
//To free the memory for struct atomic_orbital
void atomic_orbital_free(atomic_orbital * HEAD);

//To help scanning information from basis files and store it in the form of struct atomic_orbital
void basis_fscanf(FILE *,atomic_orbital *);

//set the label of a particular electron shell
void orbital_label(char *,int,int,int);
//set the angular-coefficient pair for a particular electron shell (consider d_x^2-y^2)
void orbital_angcoef_set(orbital *);

//copying the struct atomic_orbital. WARNING: you MUST allocate memory dest_HEAD before using this function
void atomic_orbital_single_cpy(atomic_orbital * dest, atomic_orbital * src);
//copying a single knot of the struct atomic_orbital. WARNING: you MUST allocate memory to dest_HEAD before using this function
void atomic_orbital_single_cpy(atomic_orbital * dest, atomic_orbital * src);
//copying the struct atomic_orbital. WARNING: you MUST allocate memory to dest_HEAD before using this function
void orbital_cpy(orbital *, orbital *);

//To count all the electron shells in this linked list
int orbital_count(orbital * HEAD);

//To synchronize the coordinates of the atom with all its electron shells
void atomic_orbital_sync_coord(atomic_orbital * atom);
//To print the name of the atom into the label of its each electron shells
void atomic_orbital_name_print(atomic_orbital * atom);


//basic functions

double double_factorial(int n);
double normalize(double alpha, int ax, int ay, int az);