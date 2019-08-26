#include "../include/RHF.h"

//obtain the quad tensor of the two-electron Coulomb integrals 
void two_electron_quad_tensor(gsl_quad_tensor * dest, orbital * HEAD, int length)
{
    int i,j,k,l;


    for(i=0;i<length;i++)
    {
        for(j=0;j<length;j++)
        {
            for(k=0;k<length;k++)
            {
                for(l=0;l<length;l++)
                {
                    gsl_quad_tensor_set(dest,i,j,k,l,two_electron_JIntegral(orbital_enquiry(HEAD,i),orbital_enquiry(HEAD,j),orbital_enquiry(HEAD,k),orbital_enquiry(HEAD,l)));
                }
            }
        }
    }
}

// calculate the Hartree-Fock energy of the system, E(HF) = \sum h_{ij} \rho_{ji} + 0.5 \sum \rho_{lk} \bar{v}_{ijkl} \rho{ji}
double HF_energy(gsl_quad_tensor * v, gsl_matrix * density_matrix, gsl_matrix * h_matrix, int length)
{
    double energy_temp = 0;

    int i,j,k,l;

    for(i=0;i<length;i++)
    {
        for(j=0;j<length;j++)
        {
            // add core hamiltonian energy
            energy_temp += gsl_matrix_get(h_matrix,i,j) * gsl_matrix_get(density_matrix,j,i);
        }
    }

    for(i=0;i<length;i++)
    {
        for(j=0;j<length;j++)
        {
            for(k=0;k<length;k++)
            {
                for(l=0;l<length;l++)
                {
                    // Coulumb integrals
                    energy_temp += gsl_quad_tensor_get(v,i,j,k,l) * gsl_matrix_get(density_matrix,j,i) * gsl_matrix_get(density_matrix,l,k);
                    // Exchange integrals
                    energy_temp -= 0.5 * gsl_quad_tensor_get(v,i,k,j,l) * gsl_matrix_get(density_matrix,j,i) * gsl_matrix_get(density_matrix,l,k);
                }
            }
        }
    }

    return 2.0 * energy_temp;
}

// calculate the attraction energy matrix (Z integrals)
double nuclear_attraction_energy_matrix_element(orbital * a, orbital * b, atomic_orbital * atom_HEAD)
{
    double result;

    atomic_orbital * atom_temp;


    result = 0;

    atom_temp = atom_HEAD;
    
    // scan through the atoms' linked list
    while(atom_temp->NEXT != NULL)
    {
        result -= atom_temp->N * orbital_ZIntegral(a,b,atom_temp->cartesian);
        atom_temp = atom_temp->NEXT;
    }

    result -= atom_temp->N * orbital_ZIntegral(a,b,atom_temp->cartesian);

    return result;
}

// calculate the single electron hamiltonian matrix (core hamiltonian matrix) element
double single_electron_hamiltonian_matrix_element(orbital * a, orbital * b, atomic_orbital * atom_HEAD)
{
    double result;

    atomic_orbital * atom_temp;

    result = orbital_kinetic_energy(a,b);

    atom_temp = atom_HEAD;
    
    while(atom_temp->NEXT != NULL)
    {
        result -= atom_temp->N * orbital_ZIntegral(a,b,atom_temp->cartesian);
        atom_temp = atom_temp->NEXT;
    }

    result -= atom_temp->N * orbital_ZIntegral(a,b,atom_temp->cartesian);

    return result;
}

// calculate the fock matrix element (F_{ij} = h_{ji} + \sum \bar{v}_{ijkl} \rho_{lk})
double fock_matrix_element(gsl_quad_tensor * v, gsl_matrix * density_matrix, gsl_matrix * h_matrix, int i, int j, int length)
{
    double fock_matrix_temp;
    int k,l;

    fock_matrix_temp = gsl_matrix_get(h_matrix,j,i);

    for(k=0;k<length;k++)
    {
        for(l=0;l<length;l++)
        {
            // Coulomb integral
            fock_matrix_temp += 2 * gsl_quad_tensor_get(v,j,i,k,l) * gsl_matrix_get(density_matrix,l,k);
            // Exchange integral 
            fock_matrix_temp -=  gsl_quad_tensor_get(v,j,k,i,l) * gsl_matrix_get(density_matrix,l,k);               
        }
    }

    return fock_matrix_temp;
}

// calculate the kinetic energy matrix 
void kinetic_energy_matrix(gsl_matrix * dest, orbital * HEAD, int length)
{
    orbital * temp1, * temp2;

    int i,j;

    i = 0;
    j = 0;

    temp1 = HEAD;
    temp2 = HEAD;
    gsl_matrix_set(dest,i,j,orbital_kinetic_energy(temp1,temp2));

    for(i=1;i<length;i++)
    {
        temp2 = temp2->NEXT;
        gsl_matrix_set(dest,i,j,orbital_kinetic_energy(temp1,temp2));
    }

    for(j=1;j<length;j++)
    {
        temp2 = HEAD;
        temp1 = temp1->NEXT;
        i = 0;
        gsl_matrix_set(dest,i,j,orbital_kinetic_energy(temp1,temp2));

        for(i=1;i<length;i++)
        {
            temp2 = temp2->NEXT;
            gsl_matrix_set(dest,i,j,orbital_kinetic_energy(temp1,temp2));
        }
    }
}

// calculate the attraction energy matrix (Z integrals)
void nuclear_attraction_energy_matrix(gsl_matrix * dest, orbital * HEAD, atomic_orbital * atom_HEAD, int length)
{
    orbital * temp1, * temp2;

    int i,j;

    i = 0;
    j = 0;

    temp1 = HEAD;
    temp2 = HEAD;

    gsl_matrix_set(dest,i,j,nuclear_attraction_energy_matrix_element(temp1,temp2,atom_HEAD));

    for(i=1;i<length;i++)
    {
        temp2 = temp2->NEXT;
        gsl_matrix_set(dest,i,j,nuclear_attraction_energy_matrix_element(temp1,temp2,atom_HEAD));
    }

    for(j=1;j<length;j++)
    {
        temp2 = HEAD;
        temp1 = temp1->NEXT;
        i = 0;
        gsl_matrix_set(dest,i,j,nuclear_attraction_energy_matrix_element(temp1,temp2,atom_HEAD));

        for(i=1;i<length;i++)
        {
            temp2 = temp2->NEXT;
            gsl_matrix_set(dest,i,j,nuclear_attraction_energy_matrix_element(temp1,temp2,atom_HEAD));
        }
    }
}

// calculate the single electron hamiltonian matrix (core hamiltonian matrix)
void core_hamiltonian_matrix(gsl_matrix * dest, orbital * HEAD, atomic_orbital * atom_HEAD, int length)
{
    orbital * temp1, * temp2;

    int i,j;

    i = 0;
    j = 0;

    temp1 = HEAD;
    temp2 = HEAD;

    gsl_matrix_set(dest,i,j,single_electron_hamiltonian_matrix_element(temp1,temp2,atom_HEAD));

    for(i=1;i<length;i++)
    {
        temp2 = temp2->NEXT;
        gsl_matrix_set(dest,i,j,single_electron_hamiltonian_matrix_element(temp1,temp2,atom_HEAD));
    }

    for(j=1;j<length;j++)
    {
        temp2 = HEAD;
        temp1 = temp1->NEXT;
        i = 0;
        gsl_matrix_set(dest,i,j,single_electron_hamiltonian_matrix_element(temp1,temp2,atom_HEAD));

        for(i=1;i<length;i++)
        {
            temp2 = temp2->NEXT;
            gsl_matrix_set(dest,i,j,single_electron_hamiltonian_matrix_element(temp1,temp2,atom_HEAD));
        }
    }    
}

//obtain the fock matrix
void fock_matrix(gsl_matrix * dest, gsl_quad_tensor * v, gsl_matrix * density_matrix, gsl_matrix * h_matrix, int length)
{
    orbital * temp1, * temp2;

    int i,j;

    for(i=0;i<length;i++)
    {
        for(j=0;j<length;j++)
        {
            gsl_matrix_set(dest,i,j,fock_matrix_element(v,density_matrix,h_matrix,i,j,length));
        }
    }
}

// perform initial guess of the coefficient matrix by performing diagonalization of core hamiltonian matrix
void initial_guess(gsl_matrix * dest, gsl_matrix * core_hamiltonian, gsl_matrix * S, int length)
{
    int i;

    i = 0;

    gsl_vector * temp_vector;

    temp_vector = gsl_vector_calloc(length);

    gsl_eigen_Lowdin_diag(core_hamiltonian,S,temp_vector,dest,length);

    gsl_vector_free(temp_vector);

}

// calculate density matrix from coefficient matrix and the number of electrons
void density_matrix(gsl_matrix * dest, gsl_matrix * coef, int el_num, int length)
{
    gsl_matrix * temp1, * temp2;

    gsl_vector * vector_temp;

    temp1 = gsl_matrix_calloc(length,length);
    temp2 = gsl_matrix_calloc(length,length);

    vector_temp = gsl_vector_calloc(length);

    int i;
    // get coefficient vectors of occupied orbitals and store into temp1, temp2
    for(i=0;i<el_num/2;i++)
    {
        gsl_matrix_get_col(vector_temp,coef,i);
        gsl_matrix_set_col(temp1,i,vector_temp);
        gsl_matrix_set_row(temp2,i,vector_temp);
    }

    // perform \sum n_a | a > < a |
    gsl_matrix_mul(temp1,temp2,dest,length,length,length);

    gsl_matrix_free(temp1);
    gsl_matrix_free(temp2);
}

//perform RHF as well as printing information
int RHF_SCF_print(double * tot_energy, gsl_vector * energy, gsl_matrix * coef, orbital * HEAD, atomic_orbital * atom_HEAD, int length, int el_num, int iteration_max, double errmax, int countmax, double alpha, int mixing_type, int SCF_INITIAL_FLAG, int SCF_FOCK_FLAG, int SCF_COEF_FLAG, int FOCK_FLAG)
{
    gsl_matrix * F, * S, * D, * h_matrix, * input_coef_temp, * output_coef_temp, * diff, * diff_temp, * mixing_temp, * debug_temp;
    
    gsl_vector * coef_vector, * difference, * vector_temp1, * vector_temp2;

    gsl_quad_tensor * v;

    // allocating memories

    // Fock Matrix
    F = gsl_matrix_calloc(length,length);
    // Overlap Matrix
    S = gsl_matrix_calloc(length,length);
    //density matrix
    D = gsl_matrix_calloc(length,length);
    // core hamiltonian matrix
    h_matrix = gsl_matrix_calloc(length,length);
    // Storing the old input coefficient matrix
    input_coef_temp = gsl_matrix_calloc(length,length);
    // Storing the old output coefficient matrix
    output_coef_temp = gsl_matrix_calloc(length,length);
    // difference between the coef and the input_coef_temp
    diff = gsl_matrix_calloc(length,length);
    // Sorting the old diff matrix
    diff_temp = gsl_matrix_calloc(length,length);
    // temporary matrix that helps mixing the coef matrix
    mixing_temp = gsl_matrix_calloc(length,length);
    // temporary matrix for debugging & showing initial information
    debug_temp = gsl_matrix_calloc(length,length);

    // vector from coef matrix
    coef_vector = gsl_vector_calloc(length);
    // norm vector of the diff matrix
    difference = gsl_vector_calloc(length);
    // temporary vector, used for mixing the coef vectors
    vector_temp1 = gsl_vector_calloc(length);
    vector_temp2 = gsl_vector_calloc(length);



    v = gsl_quad_tensor_calloc(length,length,length,length);

    int i,j,k,l,count; // iterators

    // energy variables
    double energy_temp, energy_bk;

    // parameter for mixing the previous coef vector with the new one, according to Anderson's mixing method
    double beta;
    double difference_temp;

    //Initialization
    energy_bk = 0;
    count = 0;
    beta = 0;

    printf("Nuclear repulsions: %10.6f\n",nuclei_repulsion(atom_HEAD));

    orbital_S_matrix(S,HEAD);
    core_hamiltonian_matrix(h_matrix,HEAD,atom_HEAD,length);
    initial_guess(coef,h_matrix,S,length);
    density_matrix(D,coef,el_num,length);

    if(SCF_INITIAL_FLAG==1){
        printf("Overlap Integrals:\n");
        gsl_matrix_printf(S,length,length,"%10.6f");

        kinetic_energy_matrix(debug_temp,HEAD,length);
        printf("Kinetic Energy Integrals:\n");
        gsl_matrix_printf(debug_temp,length,length,"%10.6f");

        nuclear_attraction_energy_matrix(debug_temp,HEAD,atom_HEAD,length);
        printf("Nuclear Attraction Integrals:\n");
        gsl_matrix_printf(debug_temp,length,length,"%10.6f");

        printf("Core Hamiltonian Matrix:\n");
        gsl_matrix_printf(h_matrix,length,length,"%10.6f");

        printf("\nCoefficient matrix:\n");
        gsl_matrix_printf(coef,length,length,"%10.4f");

        printf("Initial Density Matrix:\n");
        gsl_matrix_printf(D,length,length,"%10.6f");
    }

    two_electron_quad_tensor(v,HEAD,length);
    fock_matrix(F,v,D,h_matrix,length);
    energy_temp = HF_energy(v,D,h_matrix,length);
    if(SCF_INITIAL_FLAG == 0)
    {
        printf("Fock Matrix: \n");
        gsl_matrix_printf(F,length,length,"%10.6f");
        printf("Energy: %lf\n", energy_temp);
    }

    printf("\n");
    printf("============================= Start SCF =============================\n\n");
    
    for(i=0;i<iteration_max;i++)
    {
        fock_matrix(F,v,D,h_matrix,length);
        gsl_eigen_Lowdin_diag(F,S,energy,coef,length);
        density_matrix(D,coef,el_num,length);     

        energy_temp = HF_energy(v,D,h_matrix,length);

        if(fabs(energy_temp - energy_bk) < errmax) count++;
        else count = 0;

        if(count >= countmax) break;

        printf("iteration = %d, energy = %lf\n",i+1,energy_temp);
        if(SCF_FOCK_FLAG==1)
        {
            printf("\nFock matrix:\n");
            gsl_matrix_printf(F,length,length,"%10.4f");
        }
        if(SCF_COEF_FLAG==1)
        {
            printf("\nCoefficient matrix:\n");
            gsl_matrix_printf(coef,length,length,"%10.4f");
        }


        gsl_matrix_memcpy(diff_temp,diff);

        gsl_matrix_memcpy(diff,coef);
        gsl_matrix_sub(diff,input_coef_temp);


        energy_bk = energy_temp;

        // Anderson's mixing
        if(i>3 && mixing_type==1)
        {
            for(j=0;j<length;j++)
            {
                //calculating parameter beta
                gsl_matrix_get_col(vector_temp1,diff,j);
                gsl_matrix_get_col(vector_temp2,diff_temp,j);

                difference_temp = gsl_vector_inner_product(vector_temp1,vector_temp1,length) - gsl_vector_inner_product(vector_temp1,vector_temp2,length);

                gsl_vector_sub(vector_temp2,vector_temp1);

                beta = difference_temp/gsl_vector_inner_product(vector_temp2,vector_temp2,length);

                // Set the output part |n_{out}>
                gsl_matrix_get_col(vector_temp1,coef,j);
                gsl_vector_scale(vector_temp1,1.0 - beta);
                gsl_matrix_get_col(vector_temp2,output_coef_temp,j);
                gsl_vector_scale(vector_temp2,beta);
                gsl_vector_add(vector_temp1,vector_temp2);

                // Start writing the next input matrix
                gsl_matrix_set_col(mixing_temp,j,vector_temp1);

                // Set the input part |n_{in>}
                gsl_matrix_get_col(vector_temp1,input_coef_temp,j);
                gsl_vector_scale(vector_temp1,1.0-beta);
                gsl_matrix_get_col(vector_temp2,diff_temp,j);
                gsl_vector_scale(vector_temp2,beta);
                gsl_vector_sub(vector_temp1,vector_temp2);
                gsl_matrix_get_col(vector_temp2,output_coef_temp,j);
                gsl_vector_scale(vector_temp2,beta);
                gsl_vector_add(vector_temp1,vector_temp2);                

                gsl_vector_scale(vector_temp1,1.0-alpha);
                gsl_matrix_get_col(vector_temp2,mixing_temp,j);
                gsl_vector_scale(vector_temp2,alpha);
                gsl_vector_add(vector_temp1,vector_temp2);
                gsl_matrix_set_col(mixing_temp,j,vector_temp1);
            }

            // overwrite the coef matrix
            gsl_matrix_memcpy(output_coef_temp,coef);
            gsl_matrix_memcpy(coef,mixing_temp);
        }     
        gsl_matrix_memcpy(input_coef_temp,coef);   
        
    }

    printf("\n===================================================================\n\n");
    if(i==iteration_max)
        printf("WARNING: SCF not converged.\n");
    else  
        printf("SCF converged.\n");

    *tot_energy = energy_temp;
        
        if(FOCK_FLAG == 1)
        {
            printf("\nFock matrix:\n");
            gsl_matrix_printf(F,length,length,"%10.6f");            
        }

    // free memories
    gsl_matrix_free(F);
    gsl_matrix_free(S);
    gsl_matrix_free(D);
    gsl_matrix_free(h_matrix);
    gsl_matrix_free(input_coef_temp);
    gsl_matrix_free(output_coef_temp);
    gsl_matrix_free(diff);
    gsl_matrix_free(diff_temp);
    gsl_matrix_free(mixing_temp);
    gsl_matrix_free(debug_temp);

    gsl_vector_free(coef_vector);
    gsl_vector_free(difference);
    gsl_vector_free(vector_temp1);
    gsl_vector_free(vector_temp2);

    gsl_quad_tensor_free(v);  

    if(i==iteration_max)
        return 1;
    else
        return 0;
}

//calculate the nuclei repulsion energy
double nuclei_repulsion(atomic_orbital * atomlist_HEAD)
{
    atomic_orbital * temp1, * temp2;

    double result;
    double distance;

    result = 0;
    distance = 0;

    temp1 = atomlist_HEAD;
    while(temp1->NEXT != NULL)
    {
        temp2 = atomlist_HEAD;
        while(temp2->NEXT != NULL)
        {
            if(temp2 != temp1)
            {
                distance = sqrt((temp1->cartesian[0] - temp2->cartesian[0]) * (temp1->cartesian[0] - temp2->cartesian[0]) + (temp1->cartesian[1] - temp2->cartesian[1]) * (temp1->cartesian[1] - temp2->cartesian[1]) + (temp1->cartesian[2] - temp2->cartesian[2]) * (temp1->cartesian[2] - temp2->cartesian[2]));

                result += (double) temp1->N * (double) temp2->N / distance;
            }

            temp2 = temp2->NEXT;
        }
        if(temp2 != temp1)
        {
            distance = sqrt((temp1->cartesian[0] - temp2->cartesian[0]) * (temp1->cartesian[0] - temp2->cartesian[0]) + (temp1->cartesian[1] - temp2->cartesian[1]) * (temp1->cartesian[1] - temp2->cartesian[1]) + (temp1->cartesian[2] - temp2->cartesian[2]) * (temp1->cartesian[2] - temp2->cartesian[2]));

            result += (double) temp1->N * (double) temp2->N / distance;
        }

        temp1 = temp1->NEXT;
    }

    temp2 = atomlist_HEAD;
    while(temp2->NEXT != NULL)
    {
        if(temp2 != temp1)
        {
            distance = sqrt((temp1->cartesian[0] - temp2->cartesian[0]) * (temp1->cartesian[0] - temp2->cartesian[0]) + (temp1->cartesian[1] - temp2->cartesian[1]) * (temp1->cartesian[1] - temp2->cartesian[1]) + (temp1->cartesian[2] - temp2->cartesian[2]) * (temp1->cartesian[2] - temp2->cartesian[2]));

            result += (double) temp1->N * (double) temp2->N / distance;
        }

        temp2 = temp2->NEXT;
    }
    if(temp2 != temp1)
    {
        distance = sqrt((temp1->cartesian[0] - temp2->cartesian[0]) * (temp1->cartesian[0] - temp2->cartesian[0]) + (temp1->cartesian[1] - temp2->cartesian[1]) * (temp1->cartesian[1] - temp2->cartesian[1]) + (temp1->cartesian[2] - temp2->cartesian[2]) * (temp1->cartesian[2] - temp2->cartesian[2]));

    result += (double) temp1->N * (double) temp2->N / distance;
    }

    return result/2.0;
}