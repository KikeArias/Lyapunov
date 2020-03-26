#include "memoryfun.h"
int duffing(double *v, double *f,double b, double d, double z1,double e, double pp, double qq,double fi, double w);

int runge (int n, double *y, double *f, double *x, double h,double *p_Duffing)
{
int i,j;
int xaux;
double b,d,z1,e,pp,qq,fi,w; /* Para Duffing */
double savey[n+1], phi[n+1];

b=p_Duffing[0];
d=p_Duffing[1];
z1=p_Duffing[2];
e=p_Duffing[3];
pp=p_Duffing[4];
qq=p_Duffing[5];
fi=p_Duffing[6];
w=p_Duffing[7];
 

xaux=*x;

// Paso 1
duffing(y,f,b,d,z1,e,pp,qq,fi,w);

// Paso 2
for (j=1;j<=n;j++)
    {
    savey[j]=y[j];
    phi[j]=f[j];
    y[j]=savey[j]+0.5*h*f[j];
    }
duffing(y,f,b,d,z1,e,pp,qq,fi,w);

// Paso 3
for (j=1;j<=n;j++)
    {
    phi[j]=phi[j]+2.0*f[j];
    y[j]=savey[j]+0.5*h*f[j];
    }
duffing(y,f,b,d,z1,e,pp,qq,fi,w);

// Paso 4
for (j=1;j<=n;j++)
    {
    phi[j]=phi[j]+2.0*f[j];
    y[j]=savey[j]+h*f[j];
    }

duffing(y,f,b,d,z1,e,pp,qq,fi,w);

// Paso 5
 for (j=1;j<=n;j++)
     {
     y[j]=savey[j]+(phi[j]+f[j])*h/6.0;
     }
duffing(y,f,b,d,z1,e,pp,qq,fi,w);

*x=xaux+0.5*h;

return 0;
}
