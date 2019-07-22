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
        temp2 = temp1->NEXT;
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
    strcpy(pointer->name,"NULL");
    
    //Set the cartesian coordinates to zero point
    int i;
    for(i=0;i<3;i++)    pointer->cartesian[i] = 0;
    //setting the next pointer to NULL
    pointer->NEXT = NULL;

    return pointer;    
}

//To free the memory for struct atomic_orbital
void atomic_orbital_free(atomic_orbital * HEAD)
{
    atomic_orbital* temp1, * temp2;
    temp1 = HEAD;
    while(temp1->NEXT != NULL){
        temp2 = temp1->NEXT;
        orbital_free(temp1->orbital_HEAD);
        delete temp1;
        temp1 = temp2;
    }
    delete temp1;
}

//decide the label for orbital label
void orbital_label(char * dest, int n, int l, int m)
{
    char * n_str;

    strcpy(n_str,"1");
    n_str[0] = '1' + n -1;

    strcpy(dest,"");
    strcat(dest,n_str);

    switch(l)
    {
        case 0:
        strcat(dest,"s"); 
        break;

        case 1:
            switch(m)
            {
                case 0: strcat(dest,"pz"); break;
                case -1: strcat(dest,"py"); break;
                case 1: strcat (dest,"px"); break;
            }
        break;
    }

}

void basis_fscanf(FILE * basis,atomic_orbital * HEAD)
{
    atomic_orbital *temp1;

    orbital *orbit_temp1,*orbit_temp2;

    char * str;

    
    //temporary storage of angular quantum number l
    int l_temp;

    //reference of l, so that you can see whether you need to re-count the n
    int l_ref;
    
    //considering that all atoms start with 1s, having 0 angular momentum, the l_temp can be set to -1 to make a difference. so that the counting of n_counter can be operating normally
    l_temp = -1;
    l_ref = -1;

    //temporary storage of the total number of exponents and coefficients
    int tot_temp;

    //temporary storage of angular quantum number l
    int m_temp;

    int n_counter;

    n_counter = 0;

    //temporary storage of exponents and coefficients
    double * exponents_temp;
    double * coefficients_temp;

    //loops
    int i,j;

    
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
    fscanf(basis,"%d",&l_temp);
    
    //check whether it has a change compared with the former l_temp, so to restart counting the main magnetic number n
    if(l_temp == l_ref) n_counter ++;
    else n_counter = l_temp + 1;

    //set l_ref to the current value of l_temp for next loop
    l_ref = l_temp;

    //Read total
    fscanf(basis,"%s",str);
    fscanf(basis,"%d",&tot_temp);

    //Read exponents
    fscanf(basis,"%s",str);

    exponents_temp = new double[tot_temp];

    for(i=0;i<tot_temp;i++)
        fscanf(basis,"%lf",exponents_temp + i);

    //Read coefficients
    fscanf(basis,"%s",str);

    coefficients_temp = new double[tot_temp];

    for(i=0;i<tot_temp;i++)
        fscanf(basis,"%lf",exponents_temp + i);


    //create a series of orbitals in different 'direction', namely the different magnetic quantum number
    
    for(m_temp=-l_temp;m_temp<=l_temp;m_temp++)
    {
        orbit_temp1 = orbital_calloc(tot_temp);

        orbit_temp1->L = l_temp;
        orbit_temp1->total = tot_temp;
        orbit_temp1->n = n_counter;
        orbit_temp1->m = m_temp;
    }
    //set the angular quantum number as well as the total number of the coefficients and exponents
    orbit_temp1->L = l_temp;
    orbit_temp1->total = tot_temp;
    






    fscanf(basis,"%s",str);


    // while(fscanf(basis,"%s",str)!=EOF)
    // {
    //     if(strcmp(str,""))
    // }

}