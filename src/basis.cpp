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

        case 2:
            switch(m)
            {
                case -2: strcat(dest,"dx^2-y^2");break;
                case -1: strcat(dest,"dyz");break;
                case 0: strcat(dest,"dz^2");break;
                case 1: strcat(dest,"dxz");break;
                case 2: strcat(dest,"dxy");break;
            }
    }
}

//To set the angcoef of the target according to its angular quantum number and magnetic quantum number, and should correspond with the result given by orbital_label
void orbital_angcoef_set(orbital * target)
{
    switch(target->L)
    {
        case 0: 
            target->A[0].a[0] = 0;
            target->A[0].a[1] = 0;
            target->A[0].a[2] = 0;
            target->A[0].coef = 1.0;
            target->length = 1;
            break;
        
        case 1:
            switch(target->m)
            {
                case -1:
                    target->A[0].a[0] = 0;
                    target->A[0].a[1] = 1;
                    target->A[0].a[2] = 0;
                    target->A[0].coef = 1.0;
                    target->length=1;
                    break;
                case 0:
                    target->A[0].a[0] = 0;
                    target->A[0].a[1] = 0;
                    target->A[0].a[2] = 1;
                    target->A[0].coef = 1.0;
                    target->length=1;
                    break;
                case 1:
                    target->A[0].a[0] = 1;
                    target->A[0].a[1] = 0;
                    target->A[0].a[2] = 0;
                    target->A[0].coef = 1.0;
                    target->length=1;
                    break;
            }
            break;
        
        case 2:
            switch(target->m)
            {
                case -2:
                    target->A[0].a[0] = 2;
                    target->A[0].a[1] = 0;
                    target->A[0].a[2] = 0;
                    target->A[0].coef = 1.0/sqrt(2.0);
                    target->A[1].a[0] = 0;
                    target->A[1].a[1] = 2;
                    target->A[1].a[2] = 0;
                    target->A[1].coef = - 1.0/sqrt(2.0);                    
                    target->length=2;
                    break;
                case -1:
                    target->A[0].a[0] = 0;
                    target->A[0].a[1] = 1;
                    target->A[0].a[2] = 1;
                    target->A[0].coef = 1.0;   
                    target->length=1;
                    break;                 
                case 0:
                    target->A[0].a[0] = 0;
                    target->A[0].a[1] = 0;
                    target->A[0].a[2] = 2;
                    target->A[0].coef = 1.0;   
                    target->length=1;
                    break;
                case 1:
                    target->A[0].a[0] = 1;
                    target->A[0].a[1] = 0;
                    target->A[0].a[2] = 1;
                    target->A[0].coef = 1.0;   
                    target->length=1;
                    break;
                case 2:
                    target->A[0].a[0] = 1;
                    target->A[0].a[1] = 1;
                    target->A[0].a[2] = 0;
                    target->A[0].coef = 1.0;   
                    target->length=1;
                    break;    
            }
    }
}

void basis_fscanf(FILE * basis,atomic_orbital * HEAD)
{
    atomic_orbital *temp1, *temp2;

    orbital *orbit_temp1,*orbit_temp2;

    char * str;

    
    //temporary storing of angular quantum number l
    int l_temp;

    //reference of l, so that you can see whether you need to re-count the n
    int l_ref;
    
    //considering that all atoms start with 1s, having 0 angular momentum, the l_temp can be set to -1 to make a difference. so that the counting of n_counter can be operating normally
    l_temp = -1;
    l_ref = -1;

    //temporary storing of the total number of exponents and coefficients
    int tot_temp;

    //temporary storing of angular quantum number l
    int m_temp;


    //count the n
    int n_counter;
    n_counter = 0;

    //temporary storing of exponents and coefficients
    double * exponents_temp;
    double * coefficients_temp;

    //loops
    int i,j;

    //flags
    int HEADFLAG = 0;

    int ORBITHEADFLAG = 0;

    while(fscanf(basis,"%s",str)!=EOF)
    {
        //Read N
        if(strcmp(str,"N:")==0)
        {
            // Check whether this is a HEAD, namely that this is the initial part
            if(HEADFLAG == 0)
            {
                temp1 = HEAD;
                HEADFLAG = 1;
            }
            else
            {
                temp2 = temp1;
                temp1 = atomic_orbital_calloc();
                temp2->NEXT = temp1;
            }

            //Read Name
            if(strcmp(str,"Name:")==0)
            {
                fscanf(basis,"%s",str);
                strcpy(temp1->name,str);
            }
        }

        //Start reading orbitals (electronic shells)

        //Read L
        if(strcmp(str,"L=")==0)
        {
            fscanf(basis,"%s",str);
            fscanf(basis,"%d",&l_temp);
        }
        
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
            fscanf(basis,"%lf",coefficients_temp + i);


        //create a series of orbitals in different 'direction', namely the different magnetic quantum number
        
        for(m_temp=-l_temp;m_temp<=l_temp;m_temp++)
        {
            orbit_temp1 = orbital_calloc(tot_temp);
            if(ORBITHEADFLAG == 1) 
            {
                orbit_temp2->NEXT = orbit_temp1;
                temp1->orbital_HEAD = orbit_temp1;
            }
            else ORBITHEADFLAG = 1;

            orbit_temp1->L = l_temp;
            orbit_temp1->total = tot_temp;
            orbit_temp1->n = n_counter;
            orbit_temp1->m = m_temp;

            orbital_label(orbit_temp1->label,n_counter,l_temp,m_temp);
            orbital_angcoef_set(orbit_temp1);

            for(i=0;i<tot_temp;i++)
            {
                *(orbit_temp1->exponents + i) = *(exponents_temp + i);
                *(orbit_temp1->coefficients + i) = *(coefficients_temp + i);
            }

            orbit_temp2 = orbit_temp1;
        }
    }
}