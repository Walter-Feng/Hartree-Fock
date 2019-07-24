#include "../include/integral.h"

double SIntegral(double ra[3], double rb[3], int ax, int ay, int az, int bx, int by, int bz, double alpha,double beta)
{
    if(ax<0||ay<0||az<0||bx<0||by<0||bz<0) return 0;

    else if(ax > 0) return ((alpha*ra[0]+beta*rb[0])/(alpha + beta)-ra[0])*SIntegral(ra,rb,ax-1,ay,az,bx,by,bz,alpha,beta) + (ax-1)/2.0/(alpha+beta) * SIntegral(ra,rb,ax-2,ay,az,bx,by,bz,alpha,beta) + bx/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay,az,bx-1,by,bz,alpha,beta);

    else if(ay > 0) return ((alpha*ra[1]+beta*rb[1])/(alpha + beta)-ra[1])*SIntegral(ra,rb,ax,ay-1,az,bx,by,bz,alpha,beta) + (ay-1)/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay-2,az,bx,by,bz,alpha,beta) + by/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay,az,bx,by-1,bz,alpha,beta);

    else if(az > 0) return ((alpha*ra[2]+beta*rb[2])/(alpha + beta)-ra[2])*SIntegral(ra,rb,ax,ay,az-1,bx,by,bz,alpha,beta) + (az-1)/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay,az-2,bx,by,bz,alpha,beta) + bz/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay,az,bx,by,bz-1,alpha,beta);

    else if(bx > 0) return ((alpha*ra[0]+beta*rb[0])/(alpha + beta)-rb[0])*SIntegral(ra,rb,ax-1,ay,az,bx,by,bz,alpha,beta) + ax/2.0/(alpha+beta) * SIntegral(ra,rb,ax-1,ay,az,bx,by,bz,alpha,beta) + (bx-1)/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay,az,bx-2,by,bz,alpha,beta);

    else if(by > 0) return ((alpha*ra[1]+beta*rb[1])/(alpha + beta)-rb[1])*SIntegral(ra,rb,ax,ay-1,az,bx,by,bz,alpha,beta) + ay/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay-1,az,bx,by,bz,alpha,beta) + (by-1)/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay,az,bx,by-2,bz,alpha,beta);

    else if(bz > 0) return ((alpha*ra[2]+beta*rb[2])/(alpha + beta)-rb[2])*SIntegral(ra,rb,ax,ay,az,bx,by,bz-1,alpha,beta) + az/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay,az-1,bx,by,bz,alpha,beta) + (bz-1)/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay,az,bx,by,bz-2,alpha,beta) ;

    else return sqrt(M_PI/(alpha + beta)) * M_PI/(alpha + beta) * exp(- alpha * beta /(alpha + beta) * (pow(ra[0]-rb[0],2)+pow(ra[1]-rb[1],2)+pow(ra[2]-rb[2],2)));
}