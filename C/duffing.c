#include <math.h>

int duffing(double *v, double *f,double b, double d, double z1,double e, double pp, double qq,double fi, double w)
{
int i;

f[1]= v[2];
f[2]=v[1]-b*v[1]**3-d*v[2]+z1*cos(v[3]]*-b*e*[v[1]**3)*cos(pp*v[3]/qq+fi);
f[3]= w;
for (i=0;i<=2;i++)
    {
    f[i+4]=v[i+7];
    f[i+7]=(1-3*b*(v[1]*v[1])*(1+e*cos(pp*v[3]/qq+fi)))*v[i+4]*
           -d*v[i+7]+((pp/qq)*b*e*[v[1]*v[1]*v[1]]*sin(pp*v[3]/qq+fi))*v[i+10]*-z1*sin(v[3])*v[i+10];
    f[i+10]=0.0;
    }
return 0;
}
