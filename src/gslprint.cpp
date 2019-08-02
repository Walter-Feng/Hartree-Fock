#include "../include/gslprint.h"

void gsl_matrix_fprint(FILE *stream, gsl_matrix * m, int rows, int columns, char const * format)
{
    int i;
    int j;

    fprintf(stream,"\n");
    for(i=0;i<rows;i++){
        for(j=0;j<columns;j++)
        {
            fprintf(stream,format,gsl_matrix_get(m,i,j));
        };

        fprintf(stream,"\n");
    };
    fprintf(stream,"\n");   
}
void gsl_vector_fprint(FILE * stream, gsl_vector * m,int length, char const * format)
{
    int i;

    fprintf(stream,"\n");

    for(i=0;i<length;i++)
        fprintf(stream,format,gsl_vector_get(m,i));

    fprintf(stream,"\n");    
}

void gsl_matrix_printf(gsl_matrix * m, int rows, int columns, char const * format)
{
    int i;
    int j;

    printf("\n");
    for(i=0;i<rows;i++){
        for(j=0;j<columns;j++)
        {
            printf(format,gsl_matrix_get(m,i,j));
        };

        printf("\n");
    };
    printf("\n");    
}

void gsl_vector_printf(gsl_vector * m,int length, char const * format)
{
    int i;

    printf("\n");

    for(i=0;i<length;i++)
        printf(format,gsl_vector_get(m,i));

    printf("\n");    
}

void gsl_matrix_complex_fprint(FILE * stream, gsl_matrix_complex * target, int rows, int columns, char const * format)
{
    int i,j;

    fprintf(stream,"\n");

    fprintf(stream,"Real Part:\n");

    for(i=0;i<rows;++i){
        for(j=0;j<columns;++j){
            fprintf(stream,format,GSL_REAL(gsl_matrix_complex_get(target,i,j)));
        }
        fprintf(stream,"\n");
    }

    fprintf(stream,"\n");
    fprintf(stream,"Imaginary Part:\n");
    for(i=0;i<rows;++i){
        for(j=0;j<columns;++j){
            fprintf(stream,format,GSL_IMAG(gsl_matrix_complex_get(target,i,j)));
        }
        fprintf(stream,"\n");
    }
    fprintf(stream,"\n");
}

void gsl_vector_complex_fprint(FILE * stream,gsl_vector_complex * target, int length, char const * format)
{
    int i;

    fprintf(stream,"\n");

    fprintf(stream,"Real Part:\n");

    for(i=0;i<length;++i){
            fprintf(stream,format,GSL_REAL(gsl_vector_complex_get(target,i)));
        fprintf(stream,"\n");
    }

    fprintf(stream,"\n");
    fprintf(stream,"Imaginary Part:\n");
    for(i=0;i<length;++i){
            fprintf(stream,format,GSL_IMAG(gsl_vector_complex_get(target,i)));
        fprintf(stream,"\n");
    }
        fprintf(stream,"\n");

    fprintf(stream,"\n");
}

void gsl_matrix_complex_printf(gsl_matrix_complex * target, int rows, int columns, char const * format)
{
    int i,j;

    printf("\n");

    printf("Real Part:\n");

    for(i=0;i<rows;++i){
        for(j=0;j<columns;++j){
            printf(format,GSL_REAL(gsl_matrix_complex_get(target,i,j)));
        }
        printf("\n");
    }

    printf("\n");
    printf("Imaginary Part:\n");
    for(i=0;i<rows;++i){
        for(j=0;j<columns;++j){
            printf(format,GSL_IMAG(gsl_matrix_complex_get(target,i,j)));
        }
        printf("\n");
    }
    printf("\n");
}

void gsl_vector_complex_fprint(gsl_vector_complex * target, int length, char const * format)
{
    int i;

    printf("\n");

    printf("Real Part:\n");

    for(i=0;i<length;++i){
            printf(format,GSL_REAL(gsl_vector_complex_get(target,i)));
        printf("\n");
    }

    printf("\n");
    printf("Imaginary Part:\n");
    for(i=0;i<length;++i){
            printf(format,GSL_IMAG(gsl_vector_complex_get(target,i)));
        printf("\n");
    }
        printf("\n");

    printf("\n");
}
