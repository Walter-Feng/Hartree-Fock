#include "../include/integral.h"

double Gamma(double z)
{

    const int a = 12;
    static double c_space[12];
    static double *c = NULL;
    int k;
    double accm;

    if ( c == NULL ) {
        double k1_factrl = 1.0; /* (k - 1)!*(-1)^k with 0!==1*/
        c = c_space;
        c[0] = sqrt(2.0*M_PI);
        for(k=1; k < a; k++) {
            c[k] = exp(a-k) * pow(a-k, k-0.5) / k1_factrl;
            k1_factrl *= -k;
        }
    }
    accm = c[0];
    for(k=1; k < a; k++) {
        accm += c[k] / ( z + k );
    }
    accm *= exp(-(z+a)) * pow(z+a, z+0.5); /* Gamma(z+1) */
    return accm/z;
}

double Boys(double x, int n)
{
    return 0.5* pow(x,-0.5-n) * (Gamma(0.5 + n) - gsl_sf_gamma_inc(0.5+n,x));
}

double Binomials(int n, int k)
{
    int i;
    double result = 1.0;
    if(n==k) return 1;
    if(k==0) return 1;
    if(k>n) return 0;
    for(i=0;i<k;i++)
        result *= n-i;
    for(i=1;i<=k;i++)
        result /= i;

    return result;
}

double f(int k, int a, int b, double PA, double PB)
{
    int i;
    double result;
    result = 0;
    for(i=0;i<=k;i++)
        result += Binomials(a,i) * Binomials(a,k-i) * pow(PA,a-i) * pow(PB,b-k+i);

    return result;
}

double tranformationcoefficients(int a[3], int b[3], int p[3], double PA[3], double PB[3], double xi, double AB)
{
    int i;
    double result;

    result = 1;

    for(i=0;i<3;i++)
        result *= f(p[i],a[i],b[i],PA[i],PB[i]);

    result *= exp(- xi * AB); // it might cause confusion, but AB here has already been squared, meaning AB -> |AB|^2

    return result;
}

double SIntegral(double ra[3], double rb[3], int ax, int ay, int az, int bx, int by, int bz, double alpha,double beta)
{
    if(ax<0||ay<0||az<0||bx<0||by<0||bz<0) return 0;

    else if(ax > 0) return ((alpha*ra[0]+beta*rb[0])/(alpha + beta)-ra[0])*SIntegral(ra,rb,ax-1,ay,az,bx,by,bz,alpha,beta) + (ax-1)/2.0/(alpha+beta) * SIntegral(ra,rb,ax-2,ay,az,bx,by,bz,alpha,beta) + bx/2.0/(alpha+beta) * SIntegral(ra,rb,ax-1,ay,az,bx-1,by,bz,alpha,beta);

    else if(ay > 0) return ((alpha*ra[1]+beta*rb[1])/(alpha + beta)-ra[1])*SIntegral(ra,rb,ax,ay-1,az,bx,by,bz,alpha,beta) + (ay-1)/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay-2,az,bx,by,bz,alpha,beta) + by/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay-1,az,bx,by-1,bz,alpha,beta);

    else if(az > 0) return ((alpha*ra[2]+beta*rb[2])/(alpha + beta)-ra[2])*SIntegral(ra,rb,ax,ay,az-1,bx,by,bz,alpha,beta) + (az-1)/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay,az-2,bx,by,bz,alpha,beta) + bz/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay,az-1,bx,by,bz-1,alpha,beta);

    else if(bx > 0) return ((alpha*ra[0]+beta*rb[0])/(alpha + beta)-rb[0])*SIntegral(ra,rb,ax-1,ay,az,bx-1,by,bz,alpha,beta) + ax/2.0/(alpha+beta) * SIntegral(ra,rb,ax-1,ay,az,bx-1,by,bz,alpha,beta) + (bx-1)/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay,az,bx-2,by,bz,alpha,beta);

    else if(by > 0) return ((alpha*ra[1]+beta*rb[1])/(alpha + beta)-rb[1])*SIntegral(ra,rb,ax,ay-1,az,bx,by,bz,alpha,beta) + ay/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay-1,az,bx,by-1,bz,alpha,beta) + (by-1)/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay,az,bx,by-2,bz,alpha,beta);

    else if(bz > 0) return ((alpha*ra[2]+beta*rb[2])/(alpha + beta)-rb[2])*SIntegral(ra,rb,ax,ay,az,bx,by,bz-1,alpha,beta) + az/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay,az-1,bx,by,bz-1,alpha,beta) + (bz-1)/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay,az,bx,by,bz-2,alpha,beta) ;

    else return sqrt(M_PI/(alpha + beta)) * M_PI/(alpha + beta) * exp(- alpha * beta /(alpha + beta) * (pow(ra[0]-rb[0],2)+pow(ra[1]-rb[1],2)+pow(ra[2]-rb[2],2)));
}

double JIntegral(double ra[3], double rb[3], int ax, int ay, int az, int bx, int by, int bz, double alpha,double beta, int m)
{
    double zeta = alpha + beta;
    double xi = alpha * beta / zeta;
    double P[3];
    double AB;

    int i;

    for(i=0;i<3;i++)
        P[i] = (alpha * ra[i] + beta * rb[i]) / zeta;

    AB = 0;
    for(i=0;i<3;i++)
     AB += (ra[i]-rb[i])*(ra[i]-rb[i]);

    if(ax<0||ay<0||az<0||bx<0||by<0||bz<0) return 0;

    else if(ax > 0) return ((P[0]-ra[0])*JIntegral(ra,rb,ax-1,ay,az,bx,by,bz,alpha,beta,m+1) + ax/2.0/alpha * JIntegral(ra,rb,ax-2,ay,az,bx,by,bz,alpha,beta,m)- ax /2.0/ alpha * beta/zeta * JIntegral(ra,rb,ax-2,ay,az,bx,by,bz,alpha,beta,m+1) + bx/2.0 /xi * JIntegral(ra,rb,ax-1,ay,az,bx-1,by,bz,alpha,beta,m+1));

    else if(ay > 0) return ((P[1]-ra[1])*JIntegral(ra,rb,ax,ay-1,az,bx,by,bz,alpha,beta,m+1) + ay/2.0/alpha * JIntegral(ra,rb,ax,ay-2,az,bx,by,bz,alpha,beta,m)- ay /2.0/ alpha * beta/zeta * JIntegral(ra,rb,ax,ay,az-2,bx,by,bz,alpha,beta,m+1) + by/2.0 /xi * JIntegral(ra,rb,ax,ay,az-1,bx,by,bz-1,alpha,beta,m+1));

    else if(az > 0) return ((P[2]-ra[2])*JIntegral(ra,rb,ax,ay,az-1,bx,by,bz,alpha,beta,m+1) + az/2.0/alpha * JIntegral(ra,rb,ax,ay,az-2,bx,by,bz,alpha,beta,m)- az /2.0/ alpha * beta/zeta * JIntegral(ra,rb,ax,ay,az-2,bx,by,bz,alpha,beta,m+1) + bz/2.0 /xi * JIntegral(ra,rb,ax,ay,az-1,bx,by,bz-1,alpha,beta,m+1));

    else if(bx > 0) return ((P[0]-rb[0])*JIntegral(ra,rb,ax,ay,az,bx-1,by,bz,alpha,beta,m+1) + bx/2.0/alpha * JIntegral(ra,rb,ax,ay,az,bx-2,by,bz,alpha,beta,m)- bx /2.0/ alpha * beta/zeta * JIntegral(ra,rb,ax,ay,az,bx-2,by,bz,alpha,beta,m+1) + ax/2.0 /xi * JIntegral(ra,rb,ax-1,ay,az,bx-1,by,bz,alpha,beta,m+1));

    else if(by > 0) return ((P[1]-rb[1])*JIntegral(ra,rb,ax,ay,az,bx,by-1,bz,alpha,beta,m+1) + by/2.0/alpha * JIntegral(ra,rb,ax,ay,az,bx,by-2,bz,alpha,beta,m)- by /2.0/ alpha * beta/zeta * JIntegral(ra,rb,ax,ay,az,bx,by-2,bz,alpha,beta,m+1) + ay/2.0 /xi * JIntegral(ra,rb,ax,ay-1,az,bx,by-1,bz,alpha,beta,m+1));

    else if(bz > 0) return ((P[2]-rb[2])*JIntegral(ra,rb,ax,ay,az,bx,by,bz-1,alpha,beta,m+1) + bz/2.0/alpha * JIntegral(ra,rb,ax,ay,az,bx,by,bz-2,alpha,beta,m)- bz /2.0/ alpha * beta/zeta * JIntegral(ra,rb,ax,ay,az,bx,by,bz-2,alpha,beta,m+1) + az/2.0 /xi * JIntegral(ra,rb,ax,ay,az-1,bx,by,bz-1,alpha,beta,m+1));

    else return 2.0 * pow(M_PI,2.5) / alpha / beta / sqrt(zeta) *Boys(xi * AB,m);
}

void two_electron_transform(gaussian_chain * HEAD, orbital * a, orbital * b)
{
    gaussian_chain * temp, *bk;

    temp = HEAD;

    orbital * temp1, * temp2;

    int i,j,k,l;

    int q;

    int p[3];

    double coefficients;

    double xi;
    double zeta;
    double P[3];
    double PA[3];
    double PB[3];
    double AB;

    AB = 0;
    for(i=0;i<3;i++)
        AB += (a->cartesian[i]-b->cartesian[i])*(a->cartesian[i]-b->cartesian[i]);

    for(i=0;i<a->length;i++)
    {
        for(j=0;j<a->total;j++)
        {
            for(k=0;k<b->length;k++)
            {
                for(l=0;l<b->total;l++)
                {
                    zeta = *(a->exponents + j) + *(b->exponents + l);
                    xi = *(a->exponents + j) * *(b->exponents + l) / zeta;

                    for(q=0;q<3;q++)
                    {
                        P[q] = (*(a->exponents + j) * a->cartesian[q] + *(b->exponents + l) * b->cartesian[q])/zeta;
                        PA[q] = P[q] - a->cartesian[q];
                        PB[q] = P[q] - b->cartesian[q];
                    }
                    
                    for(p[0]=0;p[0]<=(a->A[i].a[0]+b->A[k].a[0]);p[0]++)
                    {
                        for(p[1]=0;p[1]<=(a->A[i].a[1]+b->A[k].a[1]);p[1]++)
                        {
                            for(p[2]=0;p[2]<=(a->A[i].a[2]+b->A[k].a[2]);p[2]++)
                            {

                                coefficients = a->A[i].coef * *(a->coefficients + j) * b->A[i].coef * *(b->coefficients + j);
                                coefficients *= tranformationcoefficients(a->A[j].a,b->A[j].a,p,PA,PB,xi,AB);

                                for(q=0;q<3;q++)
                                {
                                    temp->R[q] = P[q];
                                    temp->a[q] = p[q];
                                    temp->exponents = zeta;
                                    temp->coefficients = coefficients;
                                }
                                
                                temp->NEXT = new gaussian_chain;
                                bk = temp;
                                temp = temp->NEXT;
                            }
                        }
                    }
                }
            }
        }
    }
    delete temp;
    bk->NEXT = NULL;
}