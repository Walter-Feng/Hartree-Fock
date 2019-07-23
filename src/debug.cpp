#include "../include/basis.h"

int main(int argc, char const *argv[])
{
    int i;

    char basis[30];

    strcpy(basis,"basis/");

    atomic_orbital * basis_HEAD;

    for(i=0;i<argc;i++)
    {
        if(strcmp(argv[i],"-basis")==0)
        strcat(basis,argv[i+1]);
        strcat(basis,".txt");
    }

    FILE * basisfile = fopen(basis,"r");

    basis_fscanf(basisfile,basis_HEAD);

    return 0;
}
