#include <stdio.h>
#include <stdlib.h>
#include "memoryfun.h"

int runge (int n, double *y, double *f, double *x, double h,int m,double *phi, double *savey);
int duffing(double *v, double *f,double b, double d, double z1,double e, double pp, double qq,double fi, double w);
int clambda(double *znorm, double *v, double *gsc, double *cum, double x, double *clambda);

int main(int argc, char **argv)
{

//     CALCULO DE EXPONENTES DE LYAPUNOV
//     INTEGRADOR: RUNGE-KUTTA 4ª ORDEN
//     (ortogonalizaci¢n: iorb veces por orbita)
//     El programa da las 4 columnas de datos eta, fi, lambda+,
//     lambda-, escritas en el fichero lyamap.dat, correpondientes al
//     oscilador de Duffing forzado y excitado parametricamente, segun
//     P RL 86, 1737 (2001)
//
//     Juan Jose Miralles Canals
//     Enrique Arias Antunez
//     12/03/2018

int N,NN;
int i,j; // indices
int iii,jjj; // indices
int norbit,nporbit,iorb,npx,npy;
double *cum,*znorm,*gsc,*clambda,*v,*f;
double *phi,*savey; // Para Runge-Kutta
double pp,qq,d,z1,e,fi;
double emin,emax;
double fimin,fimax;
double b,dospi,peri,w,x;
double h;
int nstep,icount;
int m; // Para runge
 
FILE *entrada, *salida;

if (argc!=3)
   {
   printf("Usage: ./lyap input_file output_file\n");
   return 1;
   }
    
entrada=fopen(argv[1],"r");
salida=fopen(argv[2],"w");

#ifdef DEBUG
fprintf(salida,"\nDATOS DE ENTRADA\n");
fprintf(salida,"================\n");
fscanf(entrada,"%d",&N);
fscanf(entrada,"%d",&NN);
#endif

#ifdef DEBUG
fprintf(salida,"N=%d\n",N);
fprintf(salida,"NN=%d\n",NN);
#endif

b=4;
dospi=4.*1.570796326794897;
peri=dospi/1.1;
w=1.1;

// Dimensionamiento de los vectores

cum=dvector(N+1);
znorm=dvector(N+1);
gsc=dvector(N+1);
clambda=dvector(N+1);
v=dvector(NN+1);
f=dvector(NN+1);
savey=dvector(NN+1);
phi=dvector(NN+1);


// Lectura de fichero de entrada
fscanf(entrada,"%lf %lf %lf %lf %lf %lf\n",&v[1],&v[2],&pp,&qq,&d,&z1);
fscanf(entrada,"%d %d %d\n",&norbit,&nporbit,&iorb);
fscanf(entrada,"%lf %lf %d\n",&emin,&emax,&npx);
fscanf(entrada,"%lf %lf %d\n",&fimin,&fimax,&npy);

fclose(entrada);


#ifdef DEBUG
fprintf(salida,"v[1]=%lf v[2]=%lf pp=%lf qq=%lf d=%lf z1=%lf\n",v[1],v[2],pp,qq,d,z1);
fprintf(salida,"norbid=%d nporbid=%d iorb=%d\n",norbit,nporbit,iorb);
fprintf(salida,"emin=%lf emax=%lf npx=%d\n",emin,emax,npx);
fprintf(salida,"fmin=%lf fimax=%lf npy=%d\n",fimin,fimax,npy);
#endif

h=peri/(1.*nporbit);
nstep=norbit*nporbit;
iorb=nporbit/iorb;


#ifdef DEBUG
fprintf(salida,"\n============\n");
fprintf(salida,"EN EJECUCION\n");
fprintf(salida,"============\n");
fprintf(salida,"h=%lf iorb=%d nstep=%d\n",h,iorb,nstep);
#endif

//for (iii=0;iii<=npx;iii++)
for (iii=0;iii<=0.0001;iii++)
    {
    e=emin+iii*(emax-emin)/(1.0*npx);
    #ifdef DEBUG
    fprintf(salida,"Entra iii=%i e=%lf\n",iii,e);
    #endif
//    for (jjj=0;jjj<=npy;jjj++)
    for (jjj=0;jjj<=1;jjj++)
        {
        #ifdef DEBUG
        fprintf(salida,"Entra jjj=%i\n",jjj);
        #endif
        fi=fimin+jjj*(fimax-fimin)/(1.0*npy);
        h=peri/(1.0*nporbit);
        #ifdef DEBUG
        fprintf(salida,"fi=%lf h=%lf\n",fi,h);
        #endif
        /*
        // INICIALIZAR VARIABLES
        v[3]=0;
        for (i=N+1;i<=NN;i++)
            {
            v[i]=0;
            }
        #ifdef DEBUG
        for (i=0;i<=NN;i++)
            {
            fprintf(salida,"v[%d]=%lf\n",i,v[i]);
            }
        #endif
        for (i=1;i<=N;i++)
            {
            v[(N+1)*i]=1.0;
            cum[i]=0.0;
            znorm[i]=0.0;
            gsc[i]=0.0;
            }
        #ifdef DEBUG
        for (i=1;i<=N;i++)
            {
            fprintf(salida,"v[%d]=%lf\n",(N+1)*i,v[(N+1)*i]);
            fprintf(salida,"cum[%d]=%lf\n",i,cum[i]);
            fprintf(salida,"znorm[%d]=%lf\n",i,znorm[i]);
            fprintf(salida,"gsc[%d]=%lf\n",i,gsc[i]);
            }
        #endif
        x=0;
        #ifdef DEBUF
        fprintf(salida,"x=%lf\n",x);
        #endif
        icount=0;
*/
        }
    }
fprintf(salida,"\n=============\n");
fprintf(salida,"FIN EJECUCION\n");
fprintf(salida,"=============\n");

fclose(salida);
return 0;
}
