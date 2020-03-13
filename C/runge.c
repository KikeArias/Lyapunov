int runge (int n, double *y, double *f, double *x, double h,int *m,double *phi, double *savey)
{
int j;
int xaux;

xaux=*x;

switch (*m)
       {
       case 1: // Paso 1
              break;
       case 2: // Paso 2
              for (j=1;j<=n;j++)
                  {
                  savey[j]=y[j];
                  phi[j]=f[j];
                  y[j]=savey[j]+0.5*h*f[j];
                  }
              *x=xaux+0.5*h;
              break;
       case 3: // Paso 3
              for (j=1;j<=n;j++)
                  {
                  phi[j]=phi[j]+2.0*f[j];
                  y[j]=savey[j]+0.5*h*f[j];
                  }
              break;
       case 4: // Paso 4
              for (j=1;j<=n;j++)
                  {
                  phi[j]=phi[j]+2.0*f[j];
                  y[j]=savey[j]+h*f[j];
                  }
              *x=xaux+0.5*h;
              break;
       case 5: // Paso 5
              for (j=1;j<=n;j++)
                  {
                  y[j]=savey[j]+(phi[j]+f[j])*h/6.0;
                  }
              break;
       }

return 0; 
}
