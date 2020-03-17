#include <math.h>

int f_clambda(int n, double *znorm, double *v, double *gsc, double *cum, double x, double *clambda)
{
int j,k,l;

znorm[1]=0.0;
for (j=1;j<=n;j++)
    {
    znorm[1]=znorm[1]+pow(v[n*j+1],2);
    }

znorm[1]=sqrt(znorm[1]);

for (j=1;j<=n;j++)
    {
    v[n*j+1]=v[n*j+1]/znorm[1];
    }

for (j=2;j<=n;j++)
    {
    for (k=1;k<j;k++)
        {
        gsc[k]=0.0;
        for (l=1;l<=n;l++)
            {
            gsc[k]=gsc[k]+v[n*l+j]*v[n*l+k];
            }
        }
    for (k=1;k<=n;k++)
        {
        for (l=1;l<j;l++)
            {
            v[n*k+j]=v[n*k+j]-gsc[l]*v[n*k+l];
            }
        }
    znorm[j]=0.0;
    for (k=1;k<=n;k++)
        {
        znorm[j]=znorm[j]+pow(v[n*k+j],2);
        }
    znorm[j]=sqrt(znorm[j]);

    for (k=1;k<=n;k++)
        {
        v[n*k+j]=v[n*k+j]/znorm[j];
        }
    for (k=1;k<=n;k++)
        {
        cum[k]=cum[k]+log(znorm[k]);
        clambda[k]=cum[k]/x;
        }
    }

return 0;
}
