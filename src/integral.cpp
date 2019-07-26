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

double Boys(double x,int n)
{
    return 0.5* pow(x,-0.5-n) * (Gamma(0.5 + n) - gsl_sf_gamma_inc(0.5+n,x));
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

    //template ： JIntegral(ra,rb,ax,ay,az,bx,by,bz,alpha,beta,m)
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

    else return sqrt(M_PI/(alpha + beta)) * M_PI/(alpha + beta) * exp(- alpha * beta /(alpha + beta) * (pow(ra[0]-rb[0],2)+pow(ra[1]-rb[1],2)+pow(ra[2]-rb[2],2)));
}