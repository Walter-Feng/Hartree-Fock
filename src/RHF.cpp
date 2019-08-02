#include "../include/RHF.h"
#include "../include/gslextra.h"

double fock_matrix_element(orbital * a, orbital * b, orbital * HEAD, gsl_matrix * coef, int length, int el_num)
{
    double result;

    int i,j;

    gsl_vector * coef_vector_temp;

    orbital * orbital_temp;

    coef_vector_temp = gsl_vector_calloc(length);

    result = 2.0 * orbital_kinetic_energy(a,b);

    for(i=0;i<el_num/2;i++)
    {
        gsl_matrix_get_col(coef_vector_temp,coef,i);
        orbital_temp = HEAD;
        
        j=0;

        result += 4.0 * two_electron_JIntegral(a,b,orbital_temp,orbital_temp) * gsl_vector_get(coef_vector_temp,j) * gsl_vector_get(coef_vector_temp,j);
        result -= 2.0 * two_electron_JIntegral(a,orbital_temp,b,orbital_temp) * gsl_vector_get(coef_vector_temp,j) * gsl_vector_get(coef_vector_temp,j);
        for(j=1;j<length;j++)
        {
            orbital_temp = orbital_temp->NEXT;
            result += 4.0 * two_electron_JIntegral(a,b,orbital_temp,orbital_temp) * gsl_vector_get(coef_vector_temp,j) * gsl_vector_get(coef_vector_temp,j);
            result -= 2.0 * two_electron_JIntegral(a,orbital_temp,b,orbital_temp) * gsl_vector_get(coef_vector_temp,j) * gsl_vector_get(coef_vector_temp,j);
        }
    }

    gsl_vector_free(coef_vector_temp);

    return result;
}

void fock_matrix(gsl_matrix * dest, gsl_matrix * coef, orbital * HEAD, int length, int el_num)
{
    orbital * temp1, * temp2;

    int i,j;

    i = 0;
    j = 0;

    temp1 = HEAD;
    temp2 = HEAD;
    gsl_matrix_set(dest,i,j,fock_matrix_element(temp1,temp2,HEAD,coef,length,el_num));

    for(i=1;i<length;i++)
    {
        temp2 = temp2->NEXT;
        gsl_matrix_set(dest,i,j,fock_matrix_element(temp1,temp2,HEAD,coef,length,el_num));
    }

    for(j=1;j<length;j++)
    {
        temp2 = HEAD;
        gsl_matrix_set(dest,i,j,fock_matrix_element(temp1,temp2,HEAD,coef,length,el_num));

        for(i=1;i<length;i++)
        {
            temp2 = temp2->NEXT;
            gsl_matrix_set(dest,i,j,fock_matrix_element(temp1,temp2,HEAD,coef,length,el_num));
        }
    }
}

int RHF_SCF_print(gsl_vector * energy, gsl_matrix * coef, orbital * HEAD, int length, int el_num, int iteration_max, double errmax, int countmax)
{
    gsl_matrix * F, * S;

    F = gsl_matrix_calloc(length,length);
    S = gsl_matrix_calloc(length,length);

    int i,j,count;

    double energy_temp, energy_bk;

    energy_bk = 0;
    count = 0;

    orbital_S_matrix(S,HEAD);

    printf("\n");
    printf("============================= Start SCF =============================\n\n");

    for(i=0;i<iteration_max;i++)
    {
        fock_matrix(F,coef,HEAD,length,el_num);
        gsl_eigen_Lowdin_diag(F,S,energy,coef);

        energy_temp = 0;

        for(j=0;j<el_num/2;j++)
            energy_temp += gsl_vector_get(energy,i) * 2.0;

        if(abs(energy_temp - energy_bk) < errmax) count++;
        else count = 0;

        if(count >= countmax) break;

        printf("iteration = %d, energy = %lf",i+1,energy_temp);
        energy_bk = energy_temp;
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