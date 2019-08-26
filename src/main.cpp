#include "../include/RHF.h"

#define STRING_MAX 50

int main(int argc, char const *argv[])
{

    //FILES
    FILE * basisfile = NULL;
    FILE * inputfile = NULL;
    FILE * chkfile = NULL;

    //loops
    int i,j;

    int length; // the total number of the orbitals -> the width/height of the matrices
    int el_num; // the total number of the electrons
    int el_num_counter; // counter for the loops concerning electrons

    el_num = 0;

    char basis[STRING_MAX]; // file dir for basis
    char input[STRING_MAX]; // file dir for input file
    char chk[STRING_MAX]; // file dir for chk file
    char reader[STRING_MAX]; // reader for interpretation of input file

    strcpy(basis,"basis/");
    strcpy(input,"NULL");

    orbital * orbitals /*containing all orbitals*/, * orbital_temp;
    atomic_orbital * basis_HEAD; // containing all information from the basis
    atomic_orbital * atoms; // containing all information of the all atoms
    atomic_orbital * atoms_temp, * atoms_bk; // help copying the atoms from basis
    atomic_orbital * basis_scanner;

    gsl_vector * energy;
    gsl_vector * coef_vector; // extracting vectors from coefficient matrix

    gsl_matrix * coef; //coefficient matrix

    //Global (printing) parameters
    int CHK_FLAG = 0; // whether creating a chk file
    int FOCK_FLAG = 0; // whether printing fock matrix
    int COEF_FLAG = 0; // whether printing coefficient matrix

    //CHK parameters
    double CHK_i; // loop for CHK
    double CHK_j;
    double CHK_k;
    double CHK_step = 0.5;
    double CHK_BEGIN_R[3] = {-5, -5, -5}; // The loop's starting point for CHK
    double CHK_END_R[3] = {5,5,5}; // loop's ending point for CHK
    double CHK_value; // the grid values in CHK

    //SCF parmeters
    double ERR_MAX = 1e-6; // The highest tolerance of energy difference during SCF
    int SCF_MAX = 50; // Max iteration
    int COUNTER = 3; // counter ensuring that the result is converged (by repeating the result as much as counter indicates)
    int SCF_INITIAL_FLAG = 0; // whether printing initial status
    int SCF_FOCK_FLAG = 0; // whether printing fock matrix during SCF
    int SCF_COEF_FLAG = 0; // whether printing coefficient matrix during SCF
    int MIXING = 0; // choosing mixing method, currently only Anderson's mixing method is supported
    double ALPHA = 0.04; // parameter for mixing


    double coord_temp; // helping converting coordinates into atomic unit
    double total_energy = 0; // the total energy of the system


    for(i=0;i<argc;i++)
    {

        // interpreting arguments

        // specify input file
        if(strcmp(argv[i],"-f")==0)
            strcpy(input,argv[i+1]);

        // whether printing fock matrix
        if(strcmp(argv[i],"-F")==0)
            FOCK_FLAG = 1;
            
        // whether printing coefficient matrix
        if(strcmp(argv[i],"-C")==0)
            COEF_FLAG = 1;

        // print help on input files
        if(strcmp(argv[i],"--input-help")==0)
        {
            printf("\n A template of H2O with all the parameters specified: \n\n\n");
            printf("&BASIS 6-31g\n\n");
            printf("&SCF\n");
            printf("    &ERR 1e-6\n\n");

            printf("    &MAX 50\n");
            printf("    &COUNT 3\n");
            printf("    &ALPHA 0.04\n");
            printf("    &MIXING ANDERSON\n");
            printf("    &FOCK\n");
            printf("    &COEFF\n");
            printf("&END_SCF\n\n");

            printf("&COORD\n");
            printf("    H 0.0000 0.757706 -0.50839\n");
            printf("    H 0.0000 -0.757706 -0.50839\n");
            printf("    O 0.0000 0.00000 0.127098\n");
            printf("&END_COORD\n\n");

            printf("&PRINT\n");
            printf("    &CHK\n");
            printf("        &BEGIN  -5 -5 -5\n");
            printf("        &END    5  5  5\n");
            printf("        &GRID  0.05\n");
            printf("    &END_CHK\n\n");

            printf("    &FOCK\n");
            printf("    &COEFF\n");
            printf("&END_PRINT\n\n\n");
            printf("P.S.: the input file should have .in postfix. (To be honest, any postfix would work, in current version, but it is still suggested not to do so.)\n");
            printf("\n\n");
        }

        // Print help message
        if(strcmp(argv[i],"-h")==0)
        {
            printf("\n Simple help on the commands:\n\n");
            printf("-h : print this help message\n");
            printf("-f <filename> : specify the input file \n");
            printf("-F : print fock matrix (It can also be specified in the input file.)\n");
            printf("-C : print coefficient matrix (It can also be specified in the input file.)\n");
            printf("--input-help: print a template input file \n");

            printf("P.S.: the input file should have .in postfix. (To be honest, any postfix would work, in current version, but it is still suggested not to do so.)\n");
            printf("\n\n");

            return 0;
        }
    }

    // throwing error for not having an input file
    if(strcmp(input,"NULL")==0)
    {
        printf("HF_ERROR: NO INPUT FILE!\n");

        return 1;
    }

    inputfile = fopen(input,"r");
    if(inputfile == NULL)
    {
        printf("HF_ERROR: NO SUCH AN INPUT FILE!\n");

        return 1;
    }

    // interpret input file
    while(fscanf(inputfile,"%s",reader)!=EOF)
    {
        // read basis
        if(strcmp(reader,"&BASIS")==0)
        {
            // read the name of the basis
            fscanf(inputfile,"%s",reader);
            // concatenate it into basis dir
            strcat(basis,reader);
            // add postfix (to read)
            strcat(basis,".txt");
            // open file
            basisfile = fopen(basis,"r");
            //throw error if there is no basis
            if(basisfile==NULL)
            {
                printf("HF_ERROR: BASIS NOT SUPPORTED!");
                return 2;
            }
            //allocate memory for the basis_HEAD
	        basis_HEAD = atomic_orbital_calloc();
            //read basis!
            basis_fscanf(basisfile,basis_HEAD);
            fclose(basisfile);
        }

        //read SCF parameters
        if(strcmp(reader,"&SCF")==0)
        {
            while(fscanf(inputfile,"%s",reader)!=EOF)
            {
                //maximum energy difference tolerance
                if(strcmp(reader,"&ERR")==0)
                    fscanf(inputfile,"%lf",&ERR_MAX);
                //maximum SCF iteration
                if(strcmp(reader,"&MAX")==0)
                    fscanf(inputfile,"%d",&SCF_MAX);
                // maximum count
                if(strcmp(reader,"&COUNT")==0)
                    fscanf(inputfile,"%d",&COUNTER);
                // parameter for mixing method
                if(strcmp(reader,"&ALPHA")==0)
                    fscanf(inputfile,"%lf",&ALPHA);
                // choose mixing method
                if(strcmp(reader,"&MIXING")==0)
                {
                    fscanf(inputfile,"%s",reader);
                    if(strcmp(reader,"ANDERSON")==0)
                        MIXING = 1;
                }
                // whether printing initial status
                if(strcmp(reader,"&INITIAL")==0)
                    SCF_INITIAL_FLAG = 1;
                // whether printing fock matrix in SCF
                if(strcmp(reader,"&FOCK")==0)
                    SCF_FOCK_FLAG = 1;
                //whether printing coefficient matrix in SCF
                if(strcmp(reader,"&COEFF")==0)
                    SCF_COEF_FLAG = 1;
                //getting out from the loop
                if(strcmp(reader,"&END_SCF")==0)
                    break;
            }
        }

        // read atoms' information
        if(strcmp(reader,"&COORD")==0)
        {
            // allocate -> set the head of the linked list
            atoms = atomic_orbital_calloc();

            // saving the location of the head
            atoms_bk = atoms;

            // read atoms
            while(fscanf(inputfile,"%s",reader)!=EOF)
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
                    // scanner for basis
                    basis_scanner = basis_HEAD;

                    // scan the basis
                    while(basis_scanner->NEXT!=NULL)
                    {
                        // whether the name matches
                        if(strcmp(basis_scanner->name,reader)==0)
                        {
                            atomic_orbital_single_cpy(atoms,basis_scanner);

                            break;
                        }

                        basis_scanner = basis_scanner->NEXT;
                    }

                    if(strcmp(basis_scanner->name,reader)==0)
                    {
                        // copy information from the basis
                        atomic_orbital_single_cpy(atoms,basis_scanner);
                        // add electrons
                        el_num += basis_scanner->N;
                        // copy coordinates & convert unit
                        for(i=0;i<3;i++)
                        {
                            fscanf(inputfile,"%lf",&coord_temp);
                            atoms->cartesian[i] = coord_temp / 0.529177210903;
                        }

                        // copy the coordinates of the atom to all of its orbitals
                        atomic_orbital_sync_coord(atoms);
                        // print the name of the atom to all its orbitals
                        atomic_orbital_name_print(atoms);

                        atoms->NEXT = atomic_orbital_calloc();

                        atoms_temp = atoms;

                        atoms = atoms->NEXT;                        
                    }

                    // no such an atom in the basis -> throw error
                    else
                    {
                        printf("HF_ERROR: atom %s is not defined in the basis.\n",reader);

                        return 3;
                    }
                    

                }
            }
        }

        //read print parameters
        if(strcmp(reader,"&PRINT")==0)
        {
            while(fscanf(inputfile,"%s",reader)!=EOF)
            {
                if(strcmp(reader,"&END_PRINT")==0)
                    break;
                
                // read parameter for chk file
                if(strcmp(reader,"&CHK")==0)
                {
                    CHK_FLAG = 1;
                    //copy the dir & file name from input file string
                    strcpy(chk,input);
                    //add postfix
                    strcat(chk,".chk");
                    //open chk file
                    chkfile = fopen(chk,"w");

                    while(fscanf(inputfile,"%s",reader)!=EOF)
                    {
                        if(strcmp(reader,"&END_CHK")==0)
                            break;
                        else
                        {
                            // set the beginning point of the loop
                            if(strcmp(reader,"&BEGIN")==0)
                            {
                                for(i=0;i<3;i++)
                                    fscanf(inputfile,"%lf",&CHK_BEGIN_R[i]);
                            }
                            // set the ending point of the loop
                            if(strcmp(reader,"&END")==0)
                            {
                                for(i=0;i<3;i++)
                                    fscanf(inputfile,"%lf",&CHK_END_R[i]);
                            }
                            // setting the step of the loop -> the density of the grids
                            if(strcmp(reader,"&GRID")==0)
                            {
                                fscanf(inputfile,"%lf",&CHK_step);
                            }
                        }
                    }
                }
                //whether printing fock matrix
                if(strcmp(reader,"&FOCK")==0)
                    FOCK_FLAG = 1;
                //whether printing coefficient matrix
                if(strcmp(reader,"&COEFF")==0)
                    COEF_FLAG = 1;

            }
        }
    }

    // setting the head for the orbitals 
    orbitals = orbital_calloc((atoms->orbital_HEAD)->total);

    // save the location of the head
    orbital_temp = orbitals;
    // set the scanner of the atoms' list to its head
    atoms_temp = atoms;

    // copy all the orbitals from all atoms into `orbitals`
    while(atoms_temp->NEXT != NULL)
    {
        orbital_cpy(orbital_temp,atoms_temp->orbital_HEAD);
        while(orbital_temp->NEXT != NULL)
            orbital_temp = orbital_temp->NEXT;

        atoms_temp = atoms_temp->NEXT;
        orbital_temp->NEXT = orbital_calloc((atoms_temp->orbital_HEAD)->total);

        orbital_temp = orbital_temp->NEXT;
    }
    // DO the copying for the last one atom
    orbital_cpy(orbital_temp,atoms_temp->orbital_HEAD);

    // count the orbitals
    length = orbital_count(orbitals);

    // allocate memory
    energy = gsl_vector_calloc(length);
    coef_vector = gsl_vector_calloc(length);
    coef = gsl_matrix_calloc(length,length);

    printf("Total electrons: %d\n",el_num);

    //RHF!
    RHF_SCF_print(&total_energy,energy,coef,orbitals,atoms,length,el_num,SCF_MAX,ERR_MAX,COUNTER,ALPHA,MIXING,SCF_INITIAL_FLAG,SCF_FOCK_FLAG,SCF_COEF_FLAG,FOCK_FLAG);

    printf("\n\n");
    printf("MOs:\n\n");

    el_num_counter = 0;

    // print each molecular orbitals
    for(i=0;i<length;i++)
    {
        printf("MO_NUM: %d,    MO_ENERGY = %lf , occ = ",i,gsl_vector_get(energy,i));

        if(el_num_counter < el_num)
            printf("2\n");

        else printf("0\n");

        el_num_counter += 2;
    }

    printf("\nTotal Energy: %lf",total_energy + nuclei_repulsion(atoms));

    printf("\n\n");

    //print the labels of each atomic orbital in the order how they participate in the matrices
    printf("MO_LABEL:\n");
    printf("[");
    orbital_temp = orbitals;
    while(orbital_temp->NEXT != NULL)
    {
        printf(" %s ,",orbital_temp->label);
        orbital_temp = orbital_temp->NEXT;
    }
    printf(" %s ]\n\n",orbital_temp->label);

    // print the coefficient matrix
    if(COEF_FLAG == 1)
    {
        printf("MO_COEFF:\n\n");
        gsl_matrix_printf(coef,length,length,"%12.6f");
    }

    // print the chk file
    if(CHK_FLAG == 1)
    {
        //loop over the grids
        for(CHK_i = CHK_BEGIN_R[0]; CHK_i < CHK_END_R[0]; CHK_i += CHK_step)
            for(CHK_j = CHK_BEGIN_R[1]; CHK_j < CHK_END_R[1]; CHK_j += CHK_step)
                for(CHK_k = CHK_BEGIN_R[2]; CHK_k < CHK_END_R[2]; CHK_k += CHK_step)
                {
                    // select the occupying orbitals
                    CHK_value = 0;
                    for(i=0;i<el_num/2;i++)
                    {
                        // get the vector from coefficient matrix
                        gsl_matrix_get_col(coef_vector,coef,i);
                        // loop over the atomic orbitals
                        for(j=0;j<length;j++)
                        {
                            CHK_value += orbital_get(orbital_enquiry(orbitals,j),CHK_i,CHK_j,CHK_k) * gsl_vector_get(coef_vector,i) * 2.0;
                        }
                    }
                    // write into the chk file
                    fprintf(chkfile,"%20.6f%20.6f%20.6f%20.6f\n",CHK_i,CHK_j,CHK_k,CHK_value);
                }
    }
    return 0;
}
