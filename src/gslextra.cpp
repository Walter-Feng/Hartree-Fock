#include "../include/gslextra.h"

void gsl_vector_complex_convert(gsl_vector * source, gsl_vector_complex * target, int length)
{
    int i;
    double temp;
    gsl_complex complextemp;
    GSL_SET_IMAG(&complextemp,0);

    for(i=0;i<length;++i)
    {
        temp = gsl_vector_get(source,i);
        GSL_SET_REAL(&complextemp,temp);
        gsl_vector_complex_set(target,i,complextemp);
    }
}

void gsl_matrix_complex_convert(gsl_matrix * source, gsl_matrix_complex * target, int rows, int columns)
{
    int i,j;
    double temp;
    gsl_complex complextemp;
    GSL_SET_IMAG(&complextemp,0);
    
    for(i=0;i<rows;++i)
    {
        for(j=0;j<columns;++j)
        {
            temp = gsl_matrix_get(source,i,j);
            GSL_SET_REAL(&complextemp,temp);
            gsl_matrix_complex_set(target,i,j,complextemp);
        }
    }
}    

void gsl_vector_complex_extract(gsl_vector_complex * source, gsl_vector * real, gsl_vector * imag, int length)
{
    int i;

    for(i=0;i<length;++i)
    {
        gsl_vector_set(real,i,GSL_REAL(gsl_vector_complex_get(source,i)));
        gsl_vector_set(imag,i,GSL_IMAG(gsl_vector_complex_get(source,i)));
    }
}

void gsl_matrix_complex_extract(gsl_matrix_complex * source, gsl_matrix * real,gsl_matrix * imag, int rows, int columns)
{
    int i,j;

    for(i=0;i<rows;++i){
        for(j=0;j<columns;++j){
            gsl_matrix_set(real,i,j,GSL_REAL(gsl_matrix_complex_get(source,i,j)));
            gsl_matrix_set(imag,i,j,GSL_IMAG(gsl_matrix_complex_get(source,i,j)));
        }
    }
}

void gsl_matrix_diag(gsl_matrix * target, gsl_vector * diag, int length)
{
    int i;

    //initialization
    gsl_matrix_set_zero(target);

    for (i = 0; i < length; ++i)
    {
        gsl_matrix_set(target,i,i,gsl_vector_get(diag,i));
    }
}

void gsl_matrix_complex_diag(gsl_matrix_complex * target, gsl_vector_complex * diag, int length)
{
    int i;

    //initialization
    gsl_matrix_complex_set_zero(target);

    for (i = 0; i < length; ++i)
    {
        gsl_matrix_complex_set(target,i,i,gsl_vector_complex_get(diag,i));
    }
}

void gsl_matrix_mul(gsl_matrix * A, gsl_matrix *B, gsl_matrix * Result,int Acolumn,int Arow,int Bcolumn)
{
    gsl_vector * temp1;
    gsl_vector * temp2;
    
    int i;
    int j;
    int k;

    double innerproduct;

    temp1 = gsl_vector_calloc(Acolumn);
    temp2 = gsl_vector_calloc(Acolumn);
    
    gsl_matrix_set_zero(Result);

    for(i=0;i<Arow;i++){
        gsl_matrix_get_row(temp1,A,i);
        for(j=0;j<Bcolumn;j++){
            gsl_matrix_get_col(temp2,B,j);
            innerproduct=0;
            for(k=0;k<Acolumn;k++)
                innerproduct += gsl_vector_get(temp1,k) * gsl_vector_get(temp2,k);
            gsl_matrix_set(Result,i,j,innerproduct);
        }
    }
}

void gsl_matrix_complex_mul(gsl_matrix_complex * A, gsl_matrix_complex *B, gsl_matrix_complex * Result,int Acolumn,int Arow,int Bcolumn)
{
    gsl_vector_complex * temp1;
    gsl_vector_complex * temp2;
    
    int i;
    int j;
    int k;

    gsl_complex innerproduct;

    temp1 = gsl_vector_complex_calloc(Acolumn);
    temp2 = gsl_vector_complex_calloc(Acolumn);
    
    gsl_matrix_complex_set_zero(Result);

    for(i=0;i<Arow;i++){
        gsl_matrix_complex_get_row(temp1,A,i);
        for(j=0;j<Bcolumn;j++){
            gsl_matrix_complex_get_col(temp2,B,j);

            GSL_SET_REAL(&innerproduct,0);
            GSL_SET_IMAG(&innerproduct,0);

            for(k=0;k<Acolumn;k++)
                gsl_complex_add(innerproduct,gsl_complex_mul(gsl_vector_complex_get(temp1,k), gsl_vector_complex_get(temp2,k)) );
            gsl_matrix_complex_set(Result,i,j,innerproduct);
        }
    }
}

double gsl_vector_inner_product(gsl_vector * A, gsl_vector * B,int length)
{
    double result;

    result = 0;

    for(int i=0;i<length;i++){
        result += gsl_vector_get(A,i) * gsl_vector_get(B,i);
    }

    return result;
}

gsl_complex gsl_vector_complex_product(gsl_vector_complex * A, gsl_vector_complex * B, int length)
{
    gsl_complex result;
    int i;

    GSL_SET_IMAG(&result,0);
    GSL_SET_REAL(&result,0);


    for(i=0;i<length;i++){
        result = gsl_complex_add(result,gsl_complex_mul(gsl_vector_complex_get(A,i),gsl_vector_complex_get(B,i)));
    }

    return result;
}

gsl_complex gsl_vector_complex_inner_product(gsl_vector_complex * A, gsl_vector_complex * B, int length)
{
    gsl_complex result;
    int i;

    GSL_SET_IMAG(&result,0);
    GSL_SET_REAL(&result,0);

    for(i=0;i<length;i++){
        result = gsl_complex_add(result,gsl_complex_mul(gsl_complex_conjugate(gsl_vector_complex_get(A,i)),gsl_vector_complex_get(B,i)));
    }

    return result;
}

void gsl_vector_transform(gsl_vector * vec,gsl_matrix * trf,int length)
{
    gsl_vector * temp1, * temp2;
    double ip;

    temp1 = gsl_vector_calloc(length);
    temp2 = gsl_vector_calloc(length);

    for(int i=0;i<length;i++){
        gsl_matrix_get_row(temp1,trf,i);
        ip = gsl_vector_inner_product(vec,temp1,length);
        gsl_vector_set(temp2,i,ip);
    }

    gsl_vector_memcpy(vec,temp2);

    gsl_vector_free(temp1);
    gsl_vector_free(temp2);
}

void gsl_vector_complex_transform(gsl_vector_complex * vec, gsl_matrix_complex * trf,int length)
{
    gsl_vector_complex * temp1, * temp2;
    gsl_complex ip;

    temp1 = gsl_vector_complex_calloc(length);
    temp2 = gsl_vector_complex_calloc(length);

    for(int i=0;i<length;i++){
        gsl_matrix_complex_get_row(temp1,trf,i);
        ip = gsl_vector_complex_product(vec,temp1,length);
        gsl_vector_complex_set(temp2,i,ip);
    }

    gsl_vector_complex_memcpy(vec,temp2);

    gsl_vector_complex_free(temp1);
    gsl_vector_complex_free(temp2);
}

void gsl_matrix_unitmatrix(gsl_matrix * m,int length)
{
    gsl_matrix_set_zero(m);

    for(int i=0;i<length;i++)
        gsl_matrix_set(m,i,i,1);
}

void gsl_matrix_complex_unitmatrix(gsl_matrix_complex * m,int length)
{
    gsl_complex a;

    GSL_SET_REAL(&a,1);
    GSL_SET_IMAG(&a,0);
    gsl_matrix_complex_set_zero(m);

    for(int i=0;i<length;i++)
        gsl_matrix_complex_set(m,i,i,a);
}

void gsl_vector_complex_conjugate(gsl_vector_complex * v,int length)
{
    gsl_complex temp;
    int i;

    for(i=0;i<length;i++){
        temp = gsl_complex_conjugate(gsl_vector_complex_get(v,i));
        gsl_vector_complex_set(v,i,temp);
    }
}

void gsl_matrix_complex_conjugate(gsl_matrix_complex * m, int rows, int columns)
{
    gsl_complex temp;
    int i;
    int j;

    for(i=0;i<rows;i++){
        for(j=0;j<columns;j++){
            temp = gsl_complex_conjugate(gsl_matrix_complex_get(m,i,j));
            gsl_matrix_complex_set(m,i,j,temp);
        }
    }
}

void gsl_eigen_Lowdin_diag(gsl_matrix * m, gsl_matrix * S, gsl_vector * eigen, gsl_matrix * eigenvec, int length)
{
    int i;

    gsl_matrix * M, * S_minus_half, * U, * S_cpy, * temp;

    gsl_vector * S_eigen, * S_eigen_minus_half;

    M = gsl_matrix_calloc(length,length);
    S_minus_half = gsl_matrix_calloc(length,length);
    S_cpy = gsl_matrix_calloc(length,length);
    U = gsl_matrix_calloc(length,length);
    temp = gsl_matrix_calloc(length,length);

    gsl_matrix_memcpy(S_cpy,S);

    S_eigen = gsl_vector_calloc(length);
    S_eigen_minus_half = gsl_vector_calloc(length);

    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(length);

    gsl_eigen_symmv(S_cpy,S_eigen,U,w);

    for(i=0;i<length;i++)
        gsl_matrix_set(S_minus_half,i,i,1.0/sqrt(gsl_vector_get(S_eigen,i)));

    gsl_matrix_mul(U,S_minus_half,temp,length,length,length);

    gsl_matrix_transpose(U);

    gsl_matrix_mul(temp,U,S_minus_half,length,length,length);

    gsl_matrix_mul(S_minus_half,m,temp,length,length,length);
    gsl_matrix_mul(temp,S_minus_half,M,length,length,length);

    gsl_eigen_symmv(M,eigen,temp,w);

    gsl_eigen_symmv_sort(eigen,temp,GSL_EIGEN_SORT_VAL_ASC);

    gsl_matrix_mul(S_minus_half,temp,eigenvec,length,length,length);

    gsl_matrix_free(M);
    gsl_matrix_free(S_minus_half);
    gsl_matrix_free(U);
    gsl_matrix_free(S_cpy);
    gsl_matrix_free(temp);
    
    gsl_vector_free(S_eigen);
    gsl_vector_free(S_eigen_minus_half);
}

void gsl_matrix_normalize(gsl_matrix * coef,int length, int columns)
{
    int i;

    double norm;

    gsl_vector * temp;

    temp = gsl_vector_calloc(length);

    for(i=0;i<columns;i++)
    {
        gsl_matrix_get_col(temp,coef,i);
        norm = gsl_vector_inner_product(temp,temp,length);
        gsl_vector_scale(temp,sqrt(1.0/norm));
        gsl_matrix_set_col(coef,i,temp);
    }
}