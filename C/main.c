#include <stdio.h>
#include <stdlib.h>
#include "memoryfun.h"

int runge (int n, double *y, double *f, double *x, double h,double *p_Duffing);
int f_clambda(int n, double *znorm, double *v, double *gsc, double *cum, double x, double *clambda);

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
//     Actualización: 26/03/2020

int N,NN;
int i,j; // indices
int iii,jjj; // indices
int m; // Variable utilizada para las llamadas a  Runge-Kutta
int norbit,nporbit,iorb,npx,npy;
double *cum,*znorm,*gsc,*clambda,*v,*f;
double pp,qq,d,z1,e,fi;
double emin,emax;
double fimin,fimax;
double b,dospi,peri,w,x;
double h;
double *p_Duffing; /* Vector con 8 doubles para pasar como parámetros a runge que se los pase a duffing sin 
                        tener una llamada a la función muy grande */
int nstep,icount;
 
FILE *entrada, *salida;

/*
if (argc!=3)
   {
   printf("Usage: ./lyap input_file output_file\n");
   return 1;
   }
*/
    
//entrada=fopen(argv[1],"r");
//entrada=fopen("entrada.txt","r");
//salida=fopen(argv[2],"w");
salida=stdout;

#ifdef DEBUG
//fprintf(salida,"\nDATOS DE ENTRADA\n");
//fprintf(salida,"================\n");
#endif
//fscanf(entrada,"%d",&N);
//fscanf(entrada,"%d",&NN);
N=3;
NN=12;

#ifdef DEBUG
//fprintf(salida,"N=%d\n",N);
//fprintf(salida,"NN=%d\n",NN);
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
p_Duffing=dvector(8);


// Lectura de fichero de entrada
//fscanf(entrada,"%lf %lf %lf %lf %lf %lf\n",&v[1],&v[2],&pp,&qq,&d,&z1);
//fscanf(entrada,"%d %d %d\n",&norbit,&nporbit,&iorb);
//fscanf(entrada,"%lf %lf %d\n",&emin,&emax,&npx);
//fscanf(entrada,"%lf %lf %d\n",&fimin,&fimax,&npy);
v[1]=0.01;
v[2]=-0.01;
pp=1;
qq=3;
d=0.154; 
z1=0.095;
norbit=100; 
nporbit=512; 
iorb=2;
emin=0.0; 
emax=1.0;
npx=1;
fimin=0; 
fimax=6.2832; 
npy=1;
//fclose(entrada);
p_Duffing[0]=b;
p_Duffing[1]=d;
p_Duffing[2]=z1;
p_Duffing[4]=pp;
p_Duffing[5]=qq;
p_Duffing[7]=w;

h=peri/(1.*nporbit);
nstep=norbit*nporbit;
iorb=nporbit/iorb;

#ifdef DEBUG
fprintf(salida,"\n============\n");
fprintf(salida,"EN EJECUCION\n");
fprintf(salida,"============\n");
#endif
icount=0;

for (iii=0;iii<=npx;iii++)
    {
    e=emin+iii*(emax-emin)/(1.0*npx);
//    for (jjj=0;jjj<=npy;jjj)
    for (jjj=0;jjj<=npy;)
        {
        fi=fimin+jjj*(fimax-fimin)/(1.0*npy);
        h=peri/(1.0*nporbit);
        p_Duffing[3]=e;
        p_Duffing[6]=fi;
        // INICIALIZAR VARIABLES
        v[3]=0;
        for (i=N+1;i<=NN;i++)
            {
            v[i]=0;
            }
        for (i=1;i<=N;i++)
            {
            v[(N+1)*i]=1.0;
            cum[i]=0.0;
            znorm[i]=0.0;
            gsc[i]=0.0;
            }
        x=0;

        runge (NN,v,f,&x,h,p_Duffing); /* Duffing se llamará desde runge. */
        icount=icount+1;
        while ((icount%iorb)!=0)
              {
              runge (NN,v,f,&x,h,p_Duffing); /* Duffing se llamará desde runge. */
              icount=icount+1;
              printf("WHILE icount=%i iorb=%i nstep=%i iii=%i jjj=%i\n",icount,iorb,nstep,iii,jjj);
              }

        /* Cuando salga, debemos saber porqué ha salido para tomar la acción oportuna */       
        if ((icount%iorb)==0) /* Si sale por que se cumple el mçodulo */
           {
           f_clambda(N,znorm,v,gsc,cum,x,clambda);
           printf("CLAMBDA icount=%i iorb=%i nstep=%i iii=%i jjj=%i\n",icount,iorb,nstep,iii,jjj);
           if (icount==nstep) /* Si sale por alcanzar nstep */
              {
              /* Imprimir valores y poner icount=0*/
              printf("IMPRIMIRR icount=%i iorb=%i nstep=%i iii=%i jjj=%i\n",icount,iorb,nstep,iii,jjj);
              fprintf(salida,"%lf %lf %lf %lf\n",e,fi,clambda[1],clambda[2]);
              icount=0;
              jjj=jjj+1;
              return 0;
              }
           }
        } /* End if jjj */
    } /* End if iii */
fprintf(salida,"\n=============\n");
fprintf(salida,"FIN EJECUCION\n");
fprintf(salida,"=============\n");

fclose(salida);

/* Liberar toda la memoria realizada con malloc */

return 0;
}
