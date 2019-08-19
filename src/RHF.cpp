#include "../include/RHF.h"

double nuclear_attraction_energy_matrix_element(orbital * a, orbital * b, atomic_orbital * atom_HEAD)
{
    double result;

    atomic_orbital * atom_temp;


    result = 0;

    atom_temp = atom_HEAD;
    
    while(atom_temp->NEXT != NULL)
    {
        result -= atom_temp->N * orbital_ZIntegral(a,b,atom_temp->cartesian);
        atom_temp = atom_temp->NEXT;
    }

    result -= atom_temp->N * orbital_ZIntegral(a,b,atom_temp->cartesian);

    return result;
}

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


double fock_matrix_element(orbital * a, orbital * b, orbital * HEAD, atomic_orbital * atom_HEAD, gsl_matrix * coef, int length, int el_num)
{
    double result;

    int i,j;

    gsl_vector * coef_vector_temp;

    orbital * orbital_temp;

    atomic_orbital * atom_temp;

    coef_vector_temp = gsl_vector_calloc(length);

    result = orbital_kinetic_energy(a,b);

    atom_temp = atom_HEAD;
    
    while(atom_temp->NEXT != NULL)
    {
        result -= atom_temp->N * orbital_ZIntegral(a,b,atom_temp->cartesian);
        atom_temp = atom_temp->NEXT;
    }

    result -= atom_temp->N * orbital_ZIntegral(a,b,atom_temp->cartesian);

    for(i=0;i<el_num/2;i++)
    {
        gsl_matrix_get_col(coef_vector_temp,coef,i);
        orbital_temp = HEAD;
        
        j=0;

        result += 2.0 * two_electron_JIntegral(a,b,orbital_temp,orbital_temp) * gsl_vector_get(coef_vector_temp,j) * gsl_vector_get(coef_vector_temp,j);
        result -= two_electron_JIntegral(a,orbital_temp,b,orbital_temp) * gsl_vector_get(coef_vector_temp,j) * gsl_vector_get(coef_vector_temp,j);
        for(j=1;j<length;j++)
        {
            orbital_temp = orbital_temp->NEXT;
            result += 2.0 * two_electron_JIntegral(a,b,orbital_temp,orbital_temp) * gsl_vector_get(coef_vector_temp,j) * gsl_vector_get(coef_vector_temp,j);
            result -= two_electron_JIntegral(a,orbital_temp,b,orbital_temp) * gsl_vector_get(coef_vector_temp,j) * gsl_vector_get(coef_vector_temp,j);
        }
    }

    gsl_vector_free(coef_vector_temp);

    return result;
}

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

void fock_matrix(gsl_matrix * dest, gsl_matrix * coef, orbital * HEAD, atomic_orbital * atom_HEAD, int length, int el_num)
{
    orbital * temp1, * temp2;

    int i,j;

    i = 0;
    j = 0;

    temp1 = HEAD;
    temp2 = HEAD;
    gsl_matrix_set(dest,i,j,fock_matrix_element(temp1,temp2,HEAD,atom_HEAD,coef,length,el_num));

    for(i=1;i<length;i++)
    {
        temp2 = temp2->NEXT;
        gsl_matrix_set(dest,i,j,fock_matrix_element(temp1,temp2,HEAD,atom_HEAD,coef,length,el_num));
    }

    for(j=1;j<length;j++)
    {
        temp2 = HEAD;
        temp1 = temp1->NEXT;
        i = 0;
        gsl_matrix_set(dest,i,j,fock_matrix_element(temp1,temp2,HEAD,atom_HEAD,coef,length,el_num));

        for(i=1;i<length;i++)
        {
            temp2 = temp2->NEXT;
            gsl_matrix_set(dest,i,j,fock_matrix_element(temp1,temp2,HEAD,atom_HEAD,coef,length,el_num));
        }
    }
}

void initial_guess(gsl_matrix * dest, gsl_matrix * S, orbital * HEAD, atomic_orbital * atom_HEAD, int length)
{
    int i;

    i = 0;

    gsl_matrix * S_minus_half;
    gsl_matrix * dest_temp;
    gsl_vector * energy;
    orbital * temp;

    S_minus_half = gsl_matrix_calloc(length,length);
    dest_temp = gsl_matrix_calloc(length,length);

    gsl_matrix_inverse_square_root(S_minus_half,S,length);

    energy = gsl_vector_calloc(length);

    temp = HEAD;

    gsl_matrix_set(dest,i,i,1);
    gsl_vector_set(energy,i,single_electron_hamiltonian_matrix_element(temp,temp,atom_HEAD));


    for(i=1;i<length;i++)
    {
        temp = temp->NEXT;
        gsl_matrix_set(dest,i,i,1);
        gsl_vector_set(energy,i,single_electron_hamiltonian_matrix_element(temp,temp,atom_HEAD));
    }

    gsl_eigen_symmv_sort(energy,dest,GSL_EIGEN_SORT_VAL_ASC);
    gsl_matrix_mul(S_minus_half,dest,dest_temp,length,length,length);

    gsl_matrix_memcpy(dest,dest_temp);

    gsl_matrix_free(dest_temp);
    gsl_matrix_free(S_minus_half);

    gsl_vector_free(energy);

}

int RHF_SCF_print(gsl_vector * energy, gsl_matrix * coef, orbital * HEAD, atomic_orbital * atom_HEAD, int length, int el_num, int iteration_max, double errmax, int countmax, double alpha)
{
    gsl_matrix * F, * S, * input_coef_temp, * output_coef_temp, * diff, * diff_temp, * mixing_temp, * debug_temp, * density_matrix;
    
    gsl_vector * coef_vector, * difference, * vector_temp1, * vector_temp2;

    gsl_quad_tensor * v;

    // allocating memories

    // Fock Matrix
    F = gsl_matrix_calloc(length,length);
    // Overlap Matrix
    S = gsl_matrix_calloc(length,length);
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

    debug_temp = gsl_matrix_calloc(length,length);

    density_matrix = gsl_matrix_calloc(length,length);
    gsl_matrix_set(density_matrix,0,0,0.5);
    gsl_matrix_set(density_matrix,1,1,0.5);
    gsl_matrix_set(density_matrix,2,2,1);
    gsl_matrix_set(density_matrix,3,3,1);
    gsl_matrix_set(density_matrix,4,4,2.0/3.0);
    gsl_matrix_set(density_matrix,5,5,2.0/3.0);
    gsl_matrix_set(density_matrix,6,6,2.0/3.0);

    gsl_matrix_memcpy(input_coef_temp,coef);

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
    printf("Overlap Integrals:\n");
    gsl_matrix_printf(S,length,length,"%10.6f");

    kinetic_energy_matrix(debug_temp,HEAD,length);
    printf("Kinetic Energy Integrals:\n");
    gsl_matrix_printf(debug_temp,length,length,"%10.6f");

    nuclear_attraction_energy_matrix(debug_temp,HEAD,atom_HEAD,length);
    printf("Nuclear Attraction Integrals:\n");
    gsl_matrix_printf(debug_temp,length,length,"%10.6f");

    core_hamiltonian_matrix(debug_temp,HEAD,atom_HEAD,length);
    printf("Core Hamiltonian Matrix:\n");
    gsl_matrix_printf(debug_temp,length,length,"%10.6f");

    printf("Initial Density Matrix:\n");
    gsl_matrix_printf(density_matrix,length,length,"%10.6f");

    two_electron_quad_tensor(v,HEAD,length);
    energy_temp = 0;
    for(i=0;i<length;i++)
    {
        for(j=0;j<length;j++)
        {
            energy_temp += gsl_matrix_get(debug_temp,i,j) * gsl_matrix_get(density_matrix,j,i);
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
                    energy_temp += gsl_quad_tensor_get(v,i,j,k,l) * gsl_matrix_get(density_matrix,j,i) * gsl_matrix_get(density_matrix,l,k);
                    energy_temp -= 0.5 * gsl_quad_tensor_get(v,i,k,j,l) * gsl_matrix_get(density_matrix,j,i) * gsl_matrix_get(density_matrix,l,k);
                }
            }
        }
    }

    double fock_matrix_temp;

    for(i=0;i<length;i++)
    {
        for(j=0;j<length;j++)
        {
            fock_matrix_temp = gsl_matrix_get(debug_temp,i,j);

            for(k=0;k<length;k++)
            {
                for(l=0;l<length;l++)
                {
                    fock_matrix_temp += gsl_quad_tensor_get(v,i,j,k,l) * gsl_matrix_get(density_matrix,l,k);
                    fock_matrix_temp -= 0.5 * gsl_quad_tensor_get(v,i,k,j,l) * gsl_matrix_get(density_matrix,l,k);               
                }
            }

            gsl_matrix_set(F,i,j,fock_matrix_temp);
        }
    }

    printf("Fock Matrix: \n");
    gsl_matrix_printf(F,length,length,"%10.6f");

    printf("Energy: %lf\n", energy_temp);

    printf("\n");
    printf("============================= Start SCF =============================\n\n");

    // initial_guess(coef,S,HEAD,atom_HEAD,length);

    gsl_matrix_unitmatrix(coef,length);

    printf("Initial guess:\n");
    printf("\nCoefficient matrix:\n");
    gsl_matrix_printf(coef,length,length,"%10.4f");
    

    for(i=0;i<=iteration_max;i++)
    {
        fock_matrix(F,coef,HEAD,atom_HEAD,length,el_num);
        gsl_eigen_Lowdin_diag(F,S,energy,coef,length);
        // keep the sign of each vector
        for(j=0;j<length;j++)
        {
            gsl_matrix_get_col(vector_temp1,coef,j);
            gsl_matrix_get_col(vector_temp2,output_coef_temp,j);
            for(k=0;k<length;k++)
            {
                if(gsl_vector_get(vector_temp1,k)*gsl_vector_get(vector_temp2,k)<0)
                {
                    gsl_vector_scale(vector_temp1,-1.0);
                    gsl_matrix_set_col(coef,j,vector_temp1);
                    break;
                }
            }
        }

        energy_temp = 0;

        for(j=0;j<el_num/2;j++)
            energy_temp += gsl_vector_get(energy,j) * 2.0;

        if(abs(energy_temp - energy_bk) < errmax) count++;
        else count = 0;

        if(count >= countmax) break;

        printf("iteration = %d, energy = %lf\n",i+1,energy_temp);
        printf("\nFock matrix:\n");
        gsl_matrix_printf(F,length,length,"%10.4f");
        printf("\nCoefficient matrix:\n");
        gsl_matrix_printf(coef,length,length,"%10.4f");


        gsl_matrix_memcpy(diff_temp,diff);

        gsl_matrix_memcpy(diff,coef);
        gsl_matrix_sub(diff,input_coef_temp);


        energy_bk = energy_temp;

        //Anderson's mixing
        // if(i>10)
        // {
        //     for(j=0;j<length;j++)
        //     {
        //         //calculating parameter beta
        //         gsl_matrix_get_col(vector_temp1,diff,j);
        //         gsl_matrix_get_col(vector_temp2,diff_temp,j);

        //         difference_temp = gsl_vector_inner_product(vector_temp1,vector_temp1,length) - gsl_vector_inner_product(vector_temp1,vector_temp2,length);

        //         gsl_vector_sub(vector_temp2,vector_temp1);

        //         beta = difference_temp/gsl_vector_inner_product(vector_temp2,vector_temp2,length);

        //         // Set the output part |n_{out}>
        //         gsl_matrix_get_col(vector_temp1,coef,j);
        //         gsl_vector_scale(vector_temp1,1.0 - beta);
        //         gsl_matrix_get_col(vector_temp2,output_coef_temp,j);
        //         gsl_vector_scale(vector_temp2,beta);
        //         gsl_vector_add(vector_temp1,vector_temp2);

        //         // Start writing the next input matrix
        //         gsl_matrix_set_col(mixing_temp,j,vector_temp1);

        //         // Set the input part |n_{in>}
        //         gsl_matrix_get_col(vector_temp1,input_coef_temp,j);
        //         gsl_vector_scale(vector_temp1,1.0-beta);
        //         gsl_matrix_get_col(vector_temp2,diff_temp,j);
        //         gsl_vector_scale(vector_temp2,beta);
        //         gsl_vector_sub(vector_temp1,vector_temp2);
        //         gsl_matrix_get_col(vector_temp2,output_coef_temp,j);
        //         gsl_vector_scale(vector_temp2,beta);
        //         gsl_vector_add(vector_temp1,vector_temp2);                

        //         gsl_vector_scale(vector_temp1,1.0-alpha);
        //         gsl_matrix_get_col(vector_temp2,mixing_temp,j);
        //         gsl_vector_scale(vector_temp2,alpha);
        //         gsl_vector_add(vector_temp1,vector_temp2);
        //         gsl_matrix_set_col(mixing_temp,j,vector_temp1);
        //     }

            // overwrite the coef matrix
            gsl_matrix_memcpy(output_coef_temp,coef);
            // gsl_matrix_memcpy(coef,mixing_temp);
        // }     
        gsl_matrix_memcpy(input_coef_temp,coef);   
    }

    printf("\n===================================================================\n\n");
    if(i==iteration_max)
    {
        printf("WARNING: SCF not converged.\n");
        return 1;
    }
    else 
    {
        printf("SCF converged.\n");
        return 0;    
    }
}

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