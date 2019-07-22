#include "../include/basis.h"

//To allocate the memory for struct orbital
orbital* orbital_calloc(size_t n)
{
    orbital * pointer;

    //allocate memory
    pointer = new orbital;

    //Initialization
    pointer->L = 0;
    pointer->m = 0;
    pointer->n = 0;

    //Set the initial name to "NULL"
    strcpy(pointer->label,"NULL");

    //allocate the memory for the coefficients and exponents
    pointer->total = n;
    pointer->coefficients = new double[n];
    pointer->exponents = new double[n];

    //setting the next pointer to NULL
    pointer->NEXT = NULL;

    return pointer;
}

// free the memory for struct orbital
void orbital_free(orbital * HEAD)
{
    orbital* temp1, * temp2;
    temp1 = HEAD;
    while(temp1->NEXT != NULL){
        temp2 = (orbital *) temp1->NEXT;
        delete temp1->exponents;
        delete temp1->coefficients;
        delete temp1;
        temp1 = temp2;
    }
    delete temp1->exponents;
    delete temp1->coefficients;
    delete temp1;
}

//To allocate the memory for struct atomic_orbital
atomic_orbital* atomic_orbital_calloc()
{
    atomic_orbital * pointer;

    //allocate memory
    pointer = new atomic_orbital;

    //Initialization
    pointer->N = 0;
    pointer->orbital_HEAD = NULL;
    //Set the initial name to "NULL"
    pointer->name = "NULL";
    
    //Set the cartesian coordinates to zero point
    int i;
    for(i=0;i<3;i++)    pointer->cartesian[i] = 0;
    //setting the next pointer to NULL
    pointer->NEXT = NULL;

    return pointer;    
}

void atomic_orbital_free(atomic_orbital * HEAD)
{
    atomic_orbital* temp1, * temp2;
    temp1 = HEAD;
    while(temp1->NEXT != NULL){
        temp2 = (atomic_orbital *) temp1->NEXT;
        delete temp1;
        temp1 = temp2;
    }
    delete temp1;
}


void basis_fscanf(FILE * basis,atomic_orbital * HEAD)
{
    atomic_orbital *temp1;

    orbital *orbit_temp1,*orbit_temp2;

    char * str;

    

    int l;
    int tot;
    
    //Read N
    fscanf(basis,"%s",str);
    fscanf(basis,"%d",& HEAD->N);

    //Read Name
    fscanf(basis,"%s",str);
    fscanf(basis,"%s",str);
    strcpy(HEAD->name,str);

    //Start reading orbitals (electronic shells)

    //Read L
    fscanf(basis,"%s",str);
    fscanf(basis,"%d",l);

    //Read total
    fscanf(basis,"%s",str);
    fscanf(basis,"%d",tot);

    //create struct orbital
    orbit_temp1 = orbital_calloc(tot);

    orbit_temp1->L = l;
    orbit_temp1->total = tot;

    fscanf(basis,"%s",str);


    while(fscanf(basis,"%s",str)!=EOF)
    {
        if(strcmp(str,""))
    }

}