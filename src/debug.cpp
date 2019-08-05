#include "../include/RHF.h"

int main(int argc, char const *argv[])
{
    int i;

    int length;
    int el_num;
    int el_num_counter;

    el_num = 0;

    char basis[30];
    char input[30];
    char reader[30];

    strcpy(basis,"basis/");
    strcpy(input,"NULL");

    orbital * orbitals, * orbital_temp;
    atomic_orbital * basis_HEAD;
    atomic_orbital * atoms, * atoms_temp, * atoms_bk;
    atomic_orbital * basis_scanner;

    gsl_vector * energy;

    gsl_matrix * coef;


    double ERR_MAX;
    int SCF_MAX;
    int COUNTER;

    double coord_temp;
    double total_energy;

    total_energy = 0;

    

    for(i=0;i<argc;i++)
    {
        if(strcmp(argv[i],"-f")==0)
        {
            strcpy(input,argv[i+1]);
        }
    }

    if(strcmp(input,"NULL")==0)
    {
        printf("NO INPUT FILE!\n");

        return 1;
    }

    FILE * inputfile = fopen(input,"r");

    while(fscanf(inputfile,"%s",reader)!=EOF)
    {
        if(strcmp(reader,"&BASIS")==0)
        {
            fscanf(inputfile,"%s",reader);
            strcat(basis,reader);
            strcat(basis,".txt");
            FILE * basisfile = fopen(basis,"r");
            basis_HEAD = atomic_orbital_calloc();
            basis_fscanf(basisfile,basis_HEAD);
        }

        if(strcmp(reader,"&SCF"))
        {
            while(fscanf(inputfile,"%s",reader))
            {
                if(strcmp(reader,"&ERR")==0)
                    fscanf(inputfile,"%lf",&ERR_MAX);
                
                if(strcmp(reader,"&MAX")==0)
                    fscanf(inputfile,"%d",&SCF_MAX);
                
                if(strcmp(reader,"&COUNT")==0)
                    fscanf(inputfile,"%d",&COUNTER);
                if(strcmp(reader,"&END_SCF")==0)
                    break;
            }
        }

        if(strcmp(reader,"&COORD")==0)
        {
            atoms = atomic_orbital_calloc();

            atoms_bk = atoms;

            while(fscanf(inputfile,"%s",reader))
            {
                if(strcmp(reader,"&END_COORD")==0)
                {
                    atomic_orbital_free(atoms);
                    atoms = atoms_bk;
                    atoms_temp->NEXT = NULL;
                    break;
                }

                else
                {
                    basis_scanner = basis_HEAD;

                    while(basis_scanner->NEXT!=NULL)
                    {
                        if(strcmp(basis_scanner->name,reader)==0)
                        {
                            atomic_orbital_single_cpy(atoms,basis_scanner);

                            break;
                        }
                    }

                    if(strcmp(basis_scanner->name,reader)==0)
                    {
                        atomic_orbital_single_cpy(atoms,basis_scanner);

                        el_num += basis_scanner->N;
                    }

                    else
                    {
                        printf("atom %s is not defined in the basis.\n",reader);

                        return 2;
                    }
                    
                    for(i=0;i<3;i++)
                    {
                        fscanf(inputfile,"%lf",&coord_temp);
                        atoms->cartesian[i] = coord_temp / 0.529;
                    }

                    atomic_orbital_sync_coord(atoms);
                    atomic_orbital_name_print(atoms);

                    atoms->NEXT = atomic_orbital_calloc();

                    atoms_temp = atoms;

                    atoms = atoms->NEXT;

                }
            }
        }
    }

    orbitals = orbital_calloc((atoms->orbital_HEAD)->total);

    orbital_temp = orbitals;
    atoms_temp = atoms;

    while(atoms_temp->NEXT != NULL)
    {
        orbital_cpy(orbital_temp,atoms_temp->orbital_HEAD);
        while(orbital_temp->NEXT != NULL)
            orbital_temp = orbitals->NEXT;

        atoms_temp = atoms_temp->NEXT;
        orbital_temp->NEXT = orbital_calloc((atoms_temp->orbital_HEAD)->total);

        orbital_temp = orbital_temp->NEXT;
    }
    orbital_cpy(orbital_temp,atoms_temp->orbital_HEAD);


    length = orbital_count(orbitals);

    energy = gsl_vector_calloc(length);
    coef = gsl_matrix_calloc(length,length);

    for(i=0;i<length;i++)
        gsl_matrix_set(coef,i,i,1.0);
    
    RHF_SCF_print(energy,coef,orbitals,length,el_num,SCF_MAX,ERR_MAX,COUNTER);

    printf("\n\n");

    printf("Energy:\n\n");

    el_num_counter = 0;

    for(i=0;i<length;i++)
    {
        printf("MO_NUM: %d,    MO_ENERGY = %lf , occ = ",i,gsl_vector_get(energy,i));

        if(el_num_counter < el_num)
        {
            printf("2\n");
            total_energy += gsl_vector_get(energy,i) * 2.0;
        }
        else printf("0\n");

        el_num_counter += 2;
    }

    printf("\nTotal Energy: %lf",total_energy - nuclei_repulsion(atoms));

    printf("\n\n");

    printf("MO_LABEL:\n");
    printf("[");
    orbital_temp = orbitals;
    while(orbital_temp->NEXT != NULL)
    {
        printf(" %s ,",orbital_temp->label);
    }
    printf(" %s ]\n\n,",orbital_temp->label);

    printf("MO_COEFF:\n\n");
    gsl_matrix_printf(coef,length,length,"%12.6f");
    return 0;
}
