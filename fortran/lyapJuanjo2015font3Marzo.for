
      program lyapJuanjo
c      
c     CALCULO DE EXPONENTES DE LYAPUNOV
c     INTEGRADOR: RUNGE-KUTTA 4§ ORDEN  
c     (ortogonalizaci¢n: iorb veces por orbita) 
c     El programa da las 4 columnas de datos eta, fi, lambda+,
c     lambda-, escritas en el fichero lyamap.dat, correpondientes al
c     oscilador de Duffing forzado y excitado parametricamente, segun
c     P RL 86, 1737 (2001)
      parameter (N=3,NN=12)
      implicit real*8(a-h,o-z)
      integer runge,norbit,nporbit,iorb,npx,npy
      dimension cum(N),znorm(N),gsc(N),clambda(N),v(NN),f(NN)
      double precision pp,qq,d,z1,e,fi
      double precision emin,emax
      double precision fmin,fimax
c     d=0.154
c     z1=0.095
      b=4.
      dospi=4.*1.570796326794897
      peri=dospi/1.1
      w=1.1
c     Se ha sacado del programa de ricardo parameter (N=3,NN=12,X1=0.)
c     y se ha introducido parameter (N=3,NN=12), definiendo x=0.
c     Ok compila bien en Fortran Force en WXP
c     Modificado 07-Ocrubre 2007
      x=0.
      icount=0
c     Abrimos fichero lectura datos entrada
      OPEN(1,FILE='entrada_lyap.dat')
      READ(1,*,end=205)v(1),v(2),pp,qq,d,z1
      READ(1,*,end=205)norbit,nporbit,iorb
      READ(1,*,end=205)emin,emax,npx
      READ(1,*,end=205)fmin,fimax,npy
 205  CLOSE(1)
      open(5,file='lyamap.dat')
      
c     La mandanga del fichero 10, sustituye a la salida por pantalla

c     open(10,file='repetido.dat')
      
      h=peri/(1.*nporbit)
      nstep=norbit*nporbit
      iorb=nporbit/iorb
c     write(*,*)'Diagrama de bifurcaci¢n'
c     write(*,*)'emin,emax,valores entre emin y emax'
c      np->ptos entre m=0 y m=1
c      read(*,*)emin,emax,npx
c     write(*,*)
c      write(*,*)'fimin,fimax,valores entre fimin y fimax'
c     np->ptos entre m=0 y m=1
c     read(*,*)fimin,fimax,npy
c     write(*,*)'CALCULANDO'
c      open(5,file='lyamap.dat')
c      h=peri/(1.*nporbit)
c     nstep=norbit*nporbit
      Do 235 iii=0,npx
      e=emin+iii*(emax-emin)/(1.*npx)
      Do 234 jjj=0,npy 
      fi=fimin+jjj*(fimax-fimin)/(1.*npy)
      h=peri/(1.*nporbit)
c     INICIALIZAR VARIABLES
      v(3)=0.
      do 102 i=N+1,NN
       v(i)=0.
 102  continue
      do 212 i=1,N
       v((N+1)*i)=1.
       cum(i)=0.
       znorm(i)=0.
       gsc(i)=0.
 212  continue
      x=0.
c     Llamada a la funcion RUNGE-KUTTA de cuarto orden
      if (iii.eq.1) then
         stop
      endif
  11  k=runge(NN, v, f, x, h )
c     si k#1,calculo de los valores de las derivadas....
      if (k.ne.1) go to 13

C     DUFFING FORZADO Y EXCITADO PARAMETRICAMENTE
      f(1)= v(2)
      f(2)=v(1)-b*v(1)**3-d*v(2)+z1*DCOS(v(3))
     *-b*e*(v(1)**3)*DCOS(pp*v(3)/qq+fi)
      f(3)= W
      do 101 i=0,2
      f(i+4) = v(i+7)
      f(i+7)=(1-3*b*(v(1)**2)*(1+e*DCOS(pp*v(3)/qq+fi)))*v(i+4)
     *-d*v(i+7)+((pp/qq)*b*e*(v(1)**3)*DSIN(pp*v(3)/qq+fi))*v(i+10)
     *-z1*DSIN(v(3))*v(i+10)
      f(i+10)= 0.
  101 continue
      
      go to 11
  13  icount=icount+1
      if (mod(icount,iorb).eq.0) then
      znorm(1)=0.
      do 30 j=1,n
	znorm(1)=znorm(1)+v(n*j+1)**2
  30  continue
      znorm(1)=sqrt(znorm(1))
      do 40 j=1,n
       v(n*j+1)=v(n*j+1)/znorm(1)
  40  continue
      do 80 j=2,n
       do 50 k=1,(j-1)
	gsc(k)=0.
	do 50 l=1,n
	 gsc(k)=gsc(k)+v(n*l+j)*v(n*l+k)
  50    continue
       do 60 k=1,n
	do 60 l=1,(j-1)
	 v(n*k+j)=v(n*k+j)-gsc(l)*v(n*k+l)
  60   continue
       znorm(j)=0.
       do 70 k=1,n
	znorm(j)=znorm(j)+v(n*k+j)**2
  70   continue
       znorm(j)=sqrt(znorm(j))
       do 80 k=1,n
	v(n*k+j)=v(n*k+j)/znorm(j)
  80   continue
       do 90 k=1,n
c        cum(k)=cum(k)+log(znorm(k))/log(2.)
	 cum(k)=cum(k)+log(znorm(k))
	 clambda(k)=cum(k)/x
  90  continue
      else
      go to 11
      endif
      if (icount.eq.nstep) then
c     write(10,204)e,fi,clambda(1),clambda(2)
      write(*,204)e,fi,clambda(1),clambda(2)
C      write(5,204)e,fi,clambda(1),clambda(2)
      icount=0
      else
      go to 11
      endif
  234 continue
  235 continue
      close(5)
c      write(*,*)'Programa finalizado'
c     Formatos de entrada y salida.....
c 202 format(I7,f10.2,2e12.5,1e10.3)
c 202 format(4f24.7)
c 204 format(2f12.5,2e12.5)
c Cambiamos formato de salida format(2f12.5,2e12.5)
c Cambiamos formato de salida format(2f12.5,2e16.5)
c 204 format(2f12.5,2e12.5)
c 204 format(2f14.7,2f24.12)
  204 format(2f14.5,2e16.5)
      end

      function runge (n, y, f, x, h )
      implicit real*8(a-h,o-z)
      integer runge
      dimension phi(50),savey(50), y(n), f(n)
      data m/0/
      m=m+1
      go to (1,2,3,4,5),m
c     paso 1
  1   runge=1
      return
c     paso 2
  2   do 22 j=1,n
      savey(j)=y(j)
      phi(j)=f(j)
  22  y(j)=savey(j)+0.5*h*f(j)
      x=x+0.5*h
      runge=1
      return
c     paso 3
   3  do 33 j=1,n
      phi(j)=phi(j)+2.0*f(j)
  33  y(j)=savey(j)+0.5*h*f(j)
      runge=1
      return
c     paso 4
   4  do 44 j=1,n
      phi(j)=phi(j)+2.0*f(j)
  44  y(j)=savey(j)+h*f(j)
      x=x+0.5*h
      runge=1
      return
c     paso 5
   5  do 55 j=1,n
  55  y(j)=savey(j)+(phi(j)+f(j))*h/6.0
      m=0.
      runge=0.
      return
      end
