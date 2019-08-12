#include "../include/integral.h"

//Calculate the complete Gamma function
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

//Calculate the Boys function
double Boys(double x, int n)
{
    if(x==0) return (1.0/(1.0+2.0 * (double) n));
    else return 0.5* pow(x,-0.5-n) * (Gamma(0.5 + n) - gsl_sf_gamma_inc(0.5+n,x));
}

//Calculate the binomials
double Binomials(int n, int k)
{
    int i;
    double result = 1.0;
    if(n==k) return 1;
    if(k==0) return 1;
    if(k<n) return 0;
    if(k<0) return 0;
    for(i=1;i<k;i++)
        result *= (double) n-i;
    for(i=1;i<=k;i++)
        result = result /(double) i;

    return result;
}

//A special function set for calculating the transformation coefficient
double f(int k, int a, int b, double PA, double PB)
{
    int i;
    double result;
    result = 0;
    
    for(i=0;i<=k;i++)
    {
        if(a<i||b-k+i<0) continue;
        if(PA==0&&a-i==0&&PB==0&&b-k+i==0) result += Binomials(a,i) * Binomials(b,k-i);
        else if(PA==0&&a-i==0) result += Binomials(a,i) * Binomials(b,k-i)* pow(PB,b-k+i);
        else if(PB==0&&b-k+i==0) result += Binomials(a,i) * Binomials(a,k-i) * pow(PA,a-i);
        else result += Binomials(a,i) * Binomials(a,k-i) * pow(PA,a-i) * pow(PB,b-k+i);
    }

    return result;
}

//Calculate the transformation coefficient that is needed to complete the transformation from two combined gaussian function to the sum of a series of gaussian functions
double tranformation_coefficient(int a[3], int b[3], int p[3], double PA[3], double PB[3], double xi, double AB)
{
    int i;
    double result;

    result = 1;

    for(i=0;i<3;i++)
        result *= f(p[i],a[i],b[i],PA[i],PB[i]);

    result *= exp(- xi * AB); // it might cause confusion, but AB here has already been squared, meaning AB -> |AB|^2

    return result;
}

//Calculate the Overlap Integral of two gaussian function
double SIntegral(double ra[3], double rb[3], int ax, int ay, int az, int bx, int by, int bz, double alpha,double beta)
{

    // if one of the angular number goes below zero, it means it will not have contribution - the same as giving derivation to a constant
    if(ax<0||ay<0||az<0||bx<0||by<0||bz<0) return 0;

    //Provide recurrence relation
    else if(ax > 0) return ((alpha*ra[0]+beta*rb[0])/(alpha + beta)-ra[0])*SIntegral(ra,rb,ax-1,ay,az,bx,by,bz,alpha,beta) + (ax-1)/2.0/(alpha+beta) * SIntegral(ra,rb,ax-2,ay,az,bx,by,bz,alpha,beta) + bx/2.0/(alpha+beta) * SIntegral(ra,rb,ax-1,ay,az,bx-1,by,bz,alpha,beta);

    else if(ay > 0) return ((alpha*ra[1]+beta*rb[1])/(alpha + beta)-ra[1])*SIntegral(ra,rb,ax,ay-1,az,bx,by,bz,alpha,beta) + (ay-1)/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay-2,az,bx,by,bz,alpha,beta) + by/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay-1,az,bx,by-1,bz,alpha,beta);

    else if(az > 0) return ((alpha*ra[2]+beta*rb[2])/(alpha + beta)-ra[2])*SIntegral(ra,rb,ax,ay,az-1,bx,by,bz,alpha,beta) + (az-1)/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay,az-2,bx,by,bz,alpha,beta) + bz/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay,az-1,bx,by,bz-1,alpha,beta);

    else if(bx > 0) return ((alpha*ra[0]+beta*rb[0])/(alpha + beta)-rb[0])*SIntegral(ra,rb,ax,ay,az,bx-1,by,bz,alpha,beta) + ax/2.0/(alpha+beta) * SIntegral(ra,rb,ax-1,ay,az,bx-1,by,bz,alpha,beta) + (bx-1)/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay,az,bx-2,by,bz,alpha,beta);

    else if(by > 0) return ((alpha*ra[1]+beta*rb[1])/(alpha + beta)-rb[1])*SIntegral(ra,rb,ax,ay,az,bx,by-1,bz,alpha,beta) + ay/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay-1,az,bx,by-1,bz,alpha,beta) + (by-1)/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay,az,bx,by-2,bz,alpha,beta);

    else if(bz > 0) return ((alpha*ra[2]+beta*rb[2])/(alpha + beta)-rb[2])*SIntegral(ra,rb,ax,ay,az,bx,by,bz-1,alpha,beta) + az/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay,az-1,bx,by,bz-1,alpha,beta) + (bz-1)/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay,az,bx,by,bz-2,alpha,beta) ;

    //giving the starting point
    else return sqrt(M_PI/(alpha + beta)) * M_PI/(alpha + beta) * exp(- alpha * beta /(alpha + beta) * (pow(ra[0]-rb[0],2)+pow(ra[1]-rb[1],2)+pow(ra[2]-rb[2],2)));
}

//Calculate the Coulomb Integral of two gaussian function
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

    // if one of the angular number goes below zero, it means it will not have contribution - the same as giving derivation to a constant
    if(ax<0||ay<0||az<0||bx<0||by<0||bz<0) return 0;

    //Provide recurrence relation
    else if(ax > 0) return ((P[0]-ra[0])*JIntegral(ra,rb,ax-1,ay,az,bx,by,bz,alpha,beta,m+1) + (ax - 1)/2.0/alpha * JIntegral(ra,rb,ax-2,ay,az,bx,by,bz,alpha,beta,m)- (ax - 1)/2.0/ alpha * beta/zeta * JIntegral(ra,rb,ax-2,ay,az,bx,by,bz,alpha,beta,m+1) + bx/2.0 /xi * JIntegral(ra,rb,ax-1,ay,az,bx-1,by,bz,alpha,beta,m+1));

    else if(ay > 0) return ((P[1]-ra[1])*JIntegral(ra,rb,ax,ay-1,az,bx,by,bz,alpha,beta,m+1) + (ay - 1)/2.0/alpha * JIntegral(ra,rb,ax,ay-2,az,bx,by,bz,alpha,beta,m)- (ay - 1)/2.0/ alpha * beta/zeta * JIntegral(ra,rb,ax,ay-2,az,bx,by,bz,alpha,beta,m+1) + by/2.0 /xi * JIntegral(ra,rb,ax,ay-1,az,bx,by-1,bz,alpha,beta,m+1));

    else if(az > 0) return ((P[2]-ra[2])*JIntegral(ra,rb,ax,ay,az-1,bx,by,bz,alpha,beta,m+1) + (az - 1)/2.0/alpha * JIntegral(ra,rb,ax,ay,az-2,bx,by,bz,alpha,beta,m)- (az - 1)/2.0/ alpha * beta/zeta * JIntegral(ra,rb,ax,ay,az-2,bx,by,bz,alpha,beta,m+1) + bz/2.0 /xi * JIntegral(ra,rb,ax,ay,az-1,bx,by,bz-1,alpha,beta,m+1));

    else if(bx > 0) return ((P[0]-rb[0])*JIntegral(ra,rb,ax,ay,az,bx-1,by,bz,alpha,beta,m+1) + (bx - 1)/2.0/alpha * JIntegral(ra,rb,ax,ay,az,bx-2,by,bz,alpha,beta,m)- (bx - 1)/2.0/ alpha * beta/zeta * JIntegral(ra,rb,ax,ay,az,bx-2,by,bz,alpha,beta,m+1) + ax/2.0 /xi * JIntegral(ra,rb,ax-1,ay,az,bx-1,by,bz,alpha,beta,m+1));

    else if(by > 0) return ((P[1]-rb[1])*JIntegral(ra,rb,ax,ay,az,bx,by-1,bz,alpha,beta,m+1) + (by - 1)/2.0/alpha * JIntegral(ra,rb,ax,ay,az,bx,by-2,bz,alpha,beta,m)- (by - 1)/2.0/ alpha * beta/zeta * JIntegral(ra,rb,ax,ay,az,bx,by-2,bz,alpha,beta,m+1) + ay/2.0 /xi * JIntegral(ra,rb,ax,ay-1,az,bx,by-1,bz,alpha,beta,m+1));

    else if(bz > 0) return ((P[2]-rb[2])*JIntegral(ra,rb,ax,ay,az,bx,by,bz-1,alpha,beta,m+1) + (bz-1)/2.0/alpha * JIntegral(ra,rb,ax,ay,az,bx,by,bz-2,alpha,beta,m)- (bz-1)/2.0/ alpha * beta/zeta * JIntegral(ra,rb,ax,ay,az,bx,by,bz-2,alpha,beta,m+1) + az/2.0 /xi * JIntegral(ra,rb,ax,ay,az-1,bx,by,bz-1,alpha,beta,m+1));

    //Set the starting point
    else return 2.0 * pow(M_PI,2.5) / alpha / beta / sqrt(zeta) *Boys(xi * AB,m);
}

double ZIntegral(double ra[3], double rb[3], double rz[3], int ax, int ay, int az, int bx, int by, int bz, double alpha, double beta, int m)
{
    double zeta = alpha + beta;
    double xi = alpha * beta / zeta;
    double P[3];
    double AB,PC;

    int i;

    for(i=0;i<3;i++)
        P[i] = (alpha * ra[i] + beta * rb[i]) / zeta;

    AB = 0;
    for(i=0;i<3;i++)
        AB += (ra[i]-rb[i])*(ra[i]-rb[i]);

    PC = 0;
    for(i=0;i<3;i++)
        PC += (P[i] - rz[i]) * (P[i] - rz[i]);

    if(ax<0||ay<0||az<0||bx<0||by<0||bz<0) return 0;

    else if(ax > 0) return beta/(alpha+beta) * (rb[0] - ra[0]) * ZIntegral(ra,rb,rz,ax-1,ay,az,bx,by,bz,alpha,beta,m) + (rz[0] - ra[0] - beta/(alpha + beta) * (rb[0] - ra[0])) * ZIntegral(ra,rb,rz,ax-1,ay,az,bx,by,bz,alpha,beta,m+1) + (double) (ax-1) / 2.0 / (alpha + beta) * (ZIntegral(ra,rb,rz,ax-2,ay,az,bx,by,bz,alpha,beta,m) - ZIntegral(ra,rb,rz,ax-2,ay,az,bx,by,bz,alpha,beta,m+1)) + (double) bx / 2.0 / (alpha + beta) * (ZIntegral(ra,rb,rz,ax-1,ay,az,bx-1,by,bz,alpha,beta,m) - ZIntegral(ra,rb,rz,ax-1,ay,az,bx-1,by,bz,alpha,beta,m + 1));

    else if(ay > 0) return beta/(alpha+beta) * (rb[1] - ra[1]) * ZIntegral(ra,rb,rz,ax,ay-1,az,bx,by,bz,alpha,beta,m) + (rz[1] - ra[1] - beta/(alpha + beta) * (rb[1] - ra[1])) * ZIntegral(ra,rb,rz,ax,ay-1,az,bx,by,bz,alpha,beta,m+1) + (double) (ay-1) / 2.0 /(alpha + beta) * (ZIntegral(ra,rb,rz,ax,ay-2,az,bx,by,bz,alpha,beta,m) - ZIntegral(ra,rb,rz,ax,ay-2,az,bx,by,bz,alpha,beta,m+1)) + (double) by / 2.0 / (alpha + beta) * (ZIntegral(ra,rb,rz,ax,ay-1,az,bx,by-1,bz,alpha,beta,m) - ZIntegral(ra,rb,rz,ax,ay-1,az,bx,by-1,bz,alpha,beta,m + 1));

    else if(az > 0) return beta/(alpha+beta) * (rb[2] - ra[2]) * ZIntegral(ra,rb,rz,ax,ay,az-1,bx,by,bz,alpha,beta,m) + (rz[2] - ra[2] - beta/(alpha + beta) * (rb[2] - ra[2])) * ZIntegral(ra,rb,rz,ax,ay,az-1,bx,by,bz,alpha,beta,m+1) + (double) (az-1) / 2.0 /(alpha + beta) *(ZIntegral(ra,rb,rz,ax,ay,az-2,bx,by,bz,alpha,beta,m) - ZIntegral(ra,rb,rz,ax,ay,az-2,bx,by,bz,alpha,beta,m+1)) + (double) bz / 2.0 / (alpha + beta) * (ZIntegral(ra,rb,rz,ax,ay,az-1,bx,by,bz-1,alpha,beta,m) - ZIntegral(ra,rb,rz,ax,ay,az-1,bx,by,bz-1,alpha,beta,m + 1));

    else if(bx > 0) return alpha/(alpha+beta) * (ra[0] - rb[0]) * ZIntegral(ra,rb,rz,ax,ay,az,bx-1,by,bz,alpha,beta,m) + (rz[0] - rb[0] - alpha/(alpha + beta) * (ra[0] - rb[0])) * ZIntegral(ra,rb,rz,ax,ay,az,bx-1,by,bz,alpha,beta,m+1) + (double) (bx-1) / 2.0 / (alpha + beta) * (ZIntegral(ra,rb,rz,ax,ay,az,bx-2,by,bz,alpha,beta,m) - ZIntegral(ra,rb,rz,ax,ay,az,bx-2,by,bz,alpha,beta,m+1)) + (double) ax / 2.0 / (alpha + beta) * (ZIntegral(ra,rb,rz,ax-1,ay,az,bx-1,by,bz,alpha,beta,m) - ZIntegral(ra,rb,rz,ax-1,ay,az,bx-1,by,bz,alpha,beta,m + 1));

    else if(by > 0) return alpha/(alpha+beta) * (ra[1] - rb[1]) * ZIntegral(ra,rb,rz,ax,ay,az,bx,by-1,bz,alpha,beta,m) + (rz[1] - rb[1] - alpha/(alpha + beta) * (ra[1] - rb[1])) * ZIntegral(ra,rb,rz,ax,ay,az,bx,by-1,bz,alpha,beta,m+1) + (double) (by-1) / 2.0 /(alpha + beta) * (ZIntegral(ra,rb,rz,ax,ay,az,bx,by-2,bz,alpha,beta,m) - ZIntegral(ra,rb,rz,ax,ay,az,bx,by-2,bz,alpha,beta,m+1)) + (double) ay / 2.0 / (alpha + beta) * (ZIntegral(ra,rb,rz,ax,ay-1,az,bx,by-1,bz,alpha,beta,m) - ZIntegral(ra,rb,rz,ax,ay-1,az,bx,by-1,bz,alpha,beta,m + 1));

    else if(bz > 0) return alpha/(alpha+beta) * (ra[2] - rb[2]) * ZIntegral(ra,rb,rz,ax,ay,az,bx,by,bz-1,alpha,beta,m) + (rz[2] - rb[2] - alpha/(alpha + beta) * (ra[2] - rb[2])) * ZIntegral(ra,rb,rz,ax,ay,az,bx,by,bz-1,alpha,beta,m+1) + (double) (bz-1) / 2.0 /(alpha + beta) *(ZIntegral(ra,rb,rz,ax,ay,az,bx,by,bz-2,alpha,beta,m) - ZIntegral(ra,rb,rz,ax,ay,az,bx,by,bz-2,alpha,beta,m+1)) + (double) az / 2.0 / (alpha + beta) * (ZIntegral(ra,rb,rz,ax,ay,az-1,bx,by,bz-1,alpha,beta,m) - ZIntegral(ra,rb,rz,ax,ay,az-1,bx,by,bz-1,alpha,beta,m + 1));

    else return 2 * M_PI / (alpha + beta) * exp(-alpha * beta /zeta * AB) * Boys((alpha * beta) * PC,m);
}


//allocate memory for struct gaussian_chain
gaussian_chain * gaussian_chain_calloc()
{
    gaussian_chain * temp;

    int i;

    temp = new gaussian_chain;

    for(i=0;i<3;i++)
    {
        temp->R[i] = 0;
        temp->a[i] = 0;
    }
    temp->exponent = 0;
    temp->coefficient = 0;

    temp->NEXT = NULL;
    return temp;
}

//free memory for struct gaussian_chain (not suitable for a single gaussian function)
void gaussian_chain_free(gaussian_chain * HEAD)
{
    gaussian_chain * temp1, *temp2;

    temp1 = HEAD;

    while(temp1->NEXT != NULL)
    {
        temp2 = temp1;
        temp1 = temp1->NEXT;
        delete temp2;
    }
    delete temp1;
}

//transform the struct orbital to struct gaussian_chain
void single_electron_transform(gaussian_chain * HEAD, orbital * a)
{
    gaussian_chain * temp, * bk;

    temp = HEAD;

    int i,j,q;

    for(i=0;i<a->length;i++)
    {
        for(j=0;j<a->total;j++)
        {
            temp->coefficient = a->A[i].coef * *(a->coefficients + j) * normalize(*(a->exponents + j),a->A[i].a[0],a->A[i].a[1],a->A[i].a[2]);
            temp->exponent = *(a->exponents + j);
            for(q=0;q<3;q++)
            {
                temp->R[q] = a->cartesian[q];
                temp->a[q] = a->A[i].a[q];
            }
            temp->NEXT = gaussian_chain_calloc();
            bk = temp;
            temp = temp->NEXT;
        }
    }
    delete temp;
    bk->NEXT = NULL;
}

//transform the two-electron orbital to struct gaussian_chain
void two_electron_transform(gaussian_chain * HEAD, orbital * a, orbital * b)
{
    gaussian_chain * temp, *bk;

    temp = HEAD;

    orbital * temp1, * temp2;

    int i,j,k,l;

    int q;

    int p[3];

    double coefficient;

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

                                coefficient = a->A[i].coef * *(a->coefficients + j) * normalize(*(a->exponents + j),a->A[i].a[0],a->A[i].a[1],a->A[i].a[2]) * b->A[k].coef * *(b->coefficients + l) * normalize(*(b->exponents + l),b->A[k].a[0],b->A[k].a[1],b->A[k].a[2]);
                                coefficient *= tranformation_coefficient(a->A[i].a,b->A[k].a,p,PA,PB,xi,AB);

                                for(q=0;q<3;q++)
                                {
                                    temp->R[q] = P[q];
                                    temp->a[q] = p[q];
                                    temp->exponent = zeta;
                                    temp->coefficient = coefficient;
                                }
                                
                                temp->NEXT = gaussian_chain_calloc();
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

//enabling overlap integrals for gaussian_chain format
double gaussian_chain_SIntegral(gaussian_chain * a, gaussian_chain * b)
{
    return a->coefficient * b->coefficient *SIntegral(a->R,b->R,a->a[0],a->a[1],a->a[2],b->a[0],b->a[1],b->a[2],a->exponent,b->exponent);
}

//enabling Coulomb integrals for gaussian_chain format
double gaussian_chain_JIntegral(gaussian_chain * a, gaussian_chain * b)
{
    return a->coefficient * b->coefficient *JIntegral(a->R,b->R,a->a[0],a->a[1],a->a[2],b->a[0],b->a[1],b->a[2],a->exponent,b->exponent,0);
}

double gaussian_chain_ZIntegral(gaussian_chain * a, gaussian_chain * b, double rz[3])
{
    return a->coefficient * b->coefficient *ZIntegral(a->R,b->R,rz,a->a[0],a->a[1],a->a[2],b->a[0],b->a[1],b->a[2],a->exponent,b->exponent,0);
}

//enabling overlap integrals for gaussian_chain format, calculating the whole chain
double gaussian_chain_full_SIntegral(gaussian_chain * a_HEAD, gaussian_chain * b_HEAD)
{
    double result = 0;
    gaussian_chain * a_temp, * b_temp;

    a_temp = a_HEAD;
    b_temp = b_HEAD;

    while(a_temp->NEXT != NULL)
    {
        b_temp = b_HEAD;
        while(b_temp->NEXT != NULL)
        {
            result += gaussian_chain_SIntegral(a_temp,b_temp);

            b_temp = b_temp->NEXT;
        }
        result += gaussian_chain_SIntegral(a_temp,b_temp);
        a_temp = a_temp->NEXT;
    }

    b_temp = b_HEAD;
    while(b_temp->NEXT != NULL)
    {
        result += gaussian_chain_SIntegral(a_temp,b_temp);

        b_temp = b_temp->NEXT;
    }

    result += gaussian_chain_SIntegral(a_temp,b_temp);

    return result;
}

//enabling Coulomb integrals for gaussian_chain format, calculating the whole chain
double gaussian_chain_full_JIntegral(gaussian_chain * a_HEAD, gaussian_chain * b_HEAD)
{
    double result;

    result = 0;
    gaussian_chain * a_temp, * b_temp;
    
    a_temp = a_HEAD;

    while(a_temp->NEXT != NULL)
    {
        b_temp = b_HEAD;
        while(b_temp->NEXT != NULL)
        {
            result += gaussian_chain_JIntegral(a_temp,b_temp);

            b_temp = b_temp->NEXT;
        }
        result += gaussian_chain_JIntegral(a_temp,b_temp);
        a_temp = a_temp->NEXT;
    }

    b_temp = b_HEAD;
    while(b_temp->NEXT != NULL)
    {
        result += gaussian_chain_JIntegral(a_temp,b_temp);
        
        b_temp = b_temp->NEXT;
    }

    result += gaussian_chain_JIntegral(a_temp,b_temp);

    return result;
}

double gaussian_chain_full_ZIntegral(gaussian_chain * a_HEAD, gaussian_chain * b_HEAD, double rz[3])
{
    double result;

    result = 0;
    gaussian_chain * a_temp, * b_temp;
    
    a_temp = a_HEAD;

    while(a_temp->NEXT != NULL)
    {
        b_temp = b_HEAD;
        while(b_temp->NEXT != NULL)
        {
            result += gaussian_chain_ZIntegral(a_temp,b_temp,rz);

            b_temp = b_temp->NEXT;
        }
        result += gaussian_chain_ZIntegral(a_temp,b_temp,rz);
        a_temp = a_temp->NEXT;
    }

    b_temp = b_HEAD;
    while(b_temp->NEXT != NULL)
    {
        result += gaussian_chain_ZIntegral(a_temp,b_temp,rz);
        
        b_temp = b_temp->NEXT;
    }

    result += gaussian_chain_ZIntegral(a_temp,b_temp,rz);

    return result;
}

double orbital_SIntegral(orbital * a, orbital * b)
{
    double result;

    result = 0;
    gaussian_chain * a_head, * b_head, * a_temp, * b_temp;
    a_head = gaussian_chain_calloc();
    b_head = gaussian_chain_calloc();

    single_electron_transform(a_head,a);
    single_electron_transform(b_head,b);

    a_temp = a_head;

    while(a_temp->NEXT != NULL)
    {
        b_temp = b_head;
        while(b_temp->NEXT != NULL)
        {
            result += gaussian_chain_SIntegral(a_temp,b_temp);

            b_temp = b_temp->NEXT;
        }
        result += gaussian_chain_SIntegral(a_temp,b_temp);
        a_temp = a_temp->NEXT;
    }

    b_temp = b_head;
    while(b_temp->NEXT != NULL)
    {
        result += gaussian_chain_SIntegral(a_temp,b_temp);

        b_temp = b_temp->NEXT;
    }

    result += gaussian_chain_SIntegral(a_temp,b_temp);

    gaussian_chain_free(a_head);
    gaussian_chain_free(b_head);

    return result;
}

double orbital_JIntegral(orbital * a, orbital * b)
{
    double result;

    result = 0;
    gaussian_chain * a_head, * b_head, * a_temp, * b_temp;
    a_head = gaussian_chain_calloc();
    b_head = gaussian_chain_calloc();

    single_electron_transform(a_head,a);
    single_electron_transform(b_head,b);

    a_temp = a_head;

    while(a_temp->NEXT != NULL)
    {
        b_temp = b_head;
        while(b_temp->NEXT != NULL)
        {
            result += gaussian_chain_JIntegral(a_temp,b_temp);

            b_temp = b_temp->NEXT;
        }
        result += gaussian_chain_JIntegral(a_temp,b_temp);
        a_temp = a_temp->NEXT;
    }

    b_temp = b_head;
    while(b_temp->NEXT != NULL)
    {
        result += gaussian_chain_JIntegral(a_temp,b_temp);
        
        b_temp = b_temp->NEXT;
    }

    result += gaussian_chain_JIntegral(a_temp,b_temp);

    gaussian_chain_free(a_head);
    gaussian_chain_free(b_head);

    return result;
}

double orbital_ZIntegral(orbital * a, orbital * b, double rz[3])
{
    double result;

    result = 0;
    gaussian_chain * a_head, * b_head, * a_temp, * b_temp;
    a_head = gaussian_chain_calloc();
    b_head = gaussian_chain_calloc();

    single_electron_transform(a_head,a);
    single_electron_transform(b_head,b);

    a_temp = a_head;

    while(a_temp->NEXT != NULL)
    {
        b_temp = b_head;
        while(b_temp->NEXT != NULL)
        {
            result += gaussian_chain_ZIntegral(a_temp,b_temp,rz);

            b_temp = b_temp->NEXT;
        }
        result += gaussian_chain_ZIntegral(a_temp,b_temp,rz);
        a_temp = a_temp->NEXT;
    }

    b_temp = b_head;
    while(b_temp->NEXT != NULL)
    {
        result += gaussian_chain_ZIntegral(a_temp,b_temp,rz);
        
        b_temp = b_temp->NEXT;
    }

    result += gaussian_chain_ZIntegral(a_temp,b_temp,rz);

    gaussian_chain_free(a_head);
    gaussian_chain_free(b_head);

    return result;
}

double two_electron_JIntegral(orbital * a_1, orbital * b_1, orbital * c_2, orbital * d_2)
{
    double result;

    gaussian_chain * temp_1, * temp_2;

    temp_1 = gaussian_chain_calloc();
    temp_2 = gaussian_chain_calloc();

    two_electron_transform(temp_1,a_1,b_1);
    two_electron_transform(temp_2,c_2,d_2);

    result = gaussian_chain_full_JIntegral(temp_1,temp_2);
    
    gaussian_chain_free(temp_1);
    gaussian_chain_free(temp_2);

    return result;
}

void gaussian_chain_derivative(gaussian_chain * dest, gaussian_chain * src, int key)
{
    gaussian_chain * temp1, * temp2;

    int i;

    temp1 = dest;
    temp2 = src;

    while(temp2->NEXT != NULL)
    {
        if(temp2->a[key] == 0)
        {
            for(i=0;i<3;i++)
            {
                temp1->a[i] = temp2->a[i];
                temp1->R[i] = temp2->R[i];
            }
            temp1->a[key] = temp2->a[key] + 1;
            temp1->coefficient = - temp2->coefficient * 2.0 * temp2->exponent;
            temp1->exponent = temp2->exponent;
            temp1->NEXT = gaussian_chain_calloc();
            temp1 = temp1->NEXT;
            temp2 = temp2->NEXT;
        }
        else
        {
            for(i=0;i<3;i++)
            {
                temp1->a[i] = temp2->a[i];
                temp1->R[i] = temp2->R[i];
            }
            temp1->a[key] = temp2->a[key] + 1;
            temp1->coefficient = - temp2->coefficient * 2.0 * temp2->exponent;
            temp1->exponent = temp2->exponent;
            temp1->NEXT = gaussian_chain_calloc();
            temp1 = temp1->NEXT;
            for(i=0;i<3;i++)
            {
                temp1->a[i] = temp2->a[i];
                temp1->R[i] = temp2->R[i];
            }
            temp1->coefficient = (double) temp2->a[key] * temp2->coefficient;
            temp1->exponent = temp2->exponent;
            temp1->a[key] = temp2->a[key] - 1;
            temp1->NEXT = gaussian_chain_calloc();
            temp1 = temp1->NEXT;
            temp2 = temp2->NEXT;
        }

    }

    if(temp2->a[key] == 0)
    {
        for(i=0;i<3;i++)
        {
            temp1->a[i] = temp2->a[i];
            temp1->R[i] = temp2->R[i];
        }
        temp1->a[key] = temp2->a[key] + 1;
        temp1->coefficient = - temp2->coefficient * 2.0 * temp2->exponent;
        temp1->exponent = temp2->exponent;
        temp1->NEXT = gaussian_chain_calloc();
        temp1 = temp1->NEXT;
        temp2 = temp2->NEXT;
    }
    else
    {
        for(i=0;i<3;i++)
        {
            temp1->a[i] = temp2->a[i];
            temp1->R[i] = temp2->R[i];
        }
        temp1->a[key] = temp2->a[key] + 1;
        temp1->coefficient = - temp2->coefficient * 2.0 * temp2->exponent;
        temp1->exponent = temp2->exponent;
        temp1->NEXT = gaussian_chain_calloc();
        temp1 = temp1->NEXT;
        for(i=0;i<3;i++)
        {
            temp1->a[i] = temp2->a[i];
            temp1->R[i] = temp2->R[i];
        }
        temp1->coefficient = (double) temp2->a[key] * temp2->coefficient;
        temp1->exponent = temp2->exponent;
        temp1->a[key] = temp2->a[key] - 1;
    }
}

void gaussian_chain_second_derivative(gaussian_chain * dest, gaussian_chain * src, int key)
{
    gaussian_chain * temp1, * temp2, * temp3;

    temp1 = dest;
    temp2 = src;
    temp3 = gaussian_chain_calloc();
    
    gaussian_chain_derivative(temp3,temp2,key);
    gaussian_chain_derivative(temp1,temp3,key);

    gaussian_chain_free(temp3);
}

void gaussian_chain_laplacian(gaussian_chain * dest, gaussian_chain * src)
{
    int i;
    gaussian_chain * temp1, * temp2;

    temp1 = dest;
    for(i=0;i<3;i++)
    {
        gaussian_chain_second_derivative(temp1,src,i);
        while(temp1->NEXT != NULL)
            temp1 = temp1->NEXT;
        temp1->NEXT = gaussian_chain_calloc();
        temp2 = temp1;
        temp1 = temp1->NEXT;
    }
    delete temp1;
    temp2->NEXT = NULL;
}

double gaussian_chain_kinetic_energy(gaussian_chain * a_HEAD, gaussian_chain * b_HEAD)
{
    double result;
    gaussian_chain * laplacian_temp;

    laplacian_temp = gaussian_chain_calloc();

    gaussian_chain_laplacian(laplacian_temp, b_HEAD);

    result = - 0.5 * gaussian_chain_full_SIntegral(a_HEAD,laplacian_temp);

    gaussian_chain_free(laplacian_temp);

    return result;
}

double orbital_kinetic_energy(orbital * a, orbital * b)
{
    double result;

    gaussian_chain * temp1, * temp2;
    temp1 = gaussian_chain_calloc();
    temp2 = gaussian_chain_calloc();

    single_electron_transform(temp1,a);
    single_electron_transform(temp2,b);

    result = gaussian_chain_kinetic_energy(temp1,temp2);

    gaussian_chain_free(temp1);
    gaussian_chain_free(temp2);

    return result;
}

void orbital_S_matrix(gsl_matrix * dest, orbital * HEAD)
{
    orbital * temp1, * temp2;

    temp1 = HEAD;
    temp2 = HEAD;

    int i,j;

    i = 0;
    j = 0;


    while(temp1->NEXT!=NULL)
    {
        while(temp2->NEXT != NULL)
        {
            gsl_matrix_set(dest,i,j,orbital_SIntegral(temp1,temp2));
            j++;
            temp2 = temp2->NEXT;
        }
        gsl_matrix_set(dest,i,j,orbital_SIntegral(temp1,temp2));
        j = 0;
        temp2 = HEAD;

        temp1 = temp1->NEXT;
        i++;
    }

    while(temp2->NEXT != NULL)
    {
        gsl_matrix_set(dest,i,j,orbital_SIntegral(temp1,temp2));
        j++;
        temp2 = temp2->NEXT;
    }
    
    gsl_matrix_set(dest,i,j,orbital_SIntegral(temp1,temp2));
}