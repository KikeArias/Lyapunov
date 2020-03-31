
      program lyapJuanjo
c      
c     CALCULO DE EXPONENTES DE LYAPUNOV
c     INTEGRADOR: RUNGE-KUTTA 4ž ORDEN  
c     (ortogonalización: iorb veces por orbita) 
c     El programa da las 4 columnas de datos eta, fi, lambda+,
c     lambda-, escritas en el fichero lyamap.dat, correpondientes al
c     oscilador de Duffing forzado y excitado parametricamente, segun
c     P RL 86, 1737 (2001)

C     Actualizacion: Enrique Arias Antunez, 31/03/2020

      parameter (N=3,NN=12)
      implicit real*8(a-h,o-z)
      integer runge,norbit,nporbit,iorb,npx,npy
      dimension cum(N),znorm(N),gsc(N),clambda(N),v(NN),f(NN)
      double precision pp,qq,d,z1,e,fi
      double precision emin,emax
      double precision fmin,fimax
      real time_begin, time_end

      call cpu_time(time_begin)
      b=4.
      dospi=4.*1.570796326794897
      peri=dospi/1.1
      w=1.1
c     Se ha sacado del programa de ricardo parameter (N=3,NN=12,X1=0.)
c     y se ha introducido parameter (N=3,NN=12), definiendo x=0.
c     Ok compila bien en Fortran Force en WXP
c     Modificado 07-Ocrubre 2007
C      x=0.
c     Abrimos fichero lectura datos entrada
      OPEN(1,FILE='entrada_lyap.dat')
      READ(1,*,end=205)v(1),v(2),pp,qq,d,z1
      READ(1,*,end=205)norbit,nporbit,iorb
      READ(1,*,end=205)emin,emax,npx
      READ(1,*,end=205)fmin,fimax,npy
 205  CLOSE(1)
      open(5,file='lyamap.dat')
      
      h=peri/(1.*nporbit)
      nstep=norbit*nporbit
      iorb=nporbit/iorb
      icount=0
      Do iii=0,npx

         e=emin+iii*(emax-emin)/(1.*npx)

         Do jjj=0,npy
 
            fi=fimin+jjj*(fimax-fimin)/(1.*npy)
            h=peri/(1.*nporbit)

c           INICIALIZAR VARIABLES
            v(3)=0.0
            do i=N+1,NN
             v(i)=0.
            end do 

            do i=1,N
               v((N+1)*i)=1.
               cum(i)=0.
               znorm(i)=0.
               gsc(i)=0.
            end do 
            
            x=0.

c           Llamada a la funcion RUNGE-KUTTA de cuarto orden
C           La variable loop controla la permanencia en el while. Gracias a este mecanismo nos evitamos
C           el uso de los goto 
            loop=0 

            Do WHILE (loop.eq.0)

               k=runge(NN, v, f, x, h )

c             si k#1,calculo de los valores de las derivadas....
              if (k.eq.1) then

C                DUFFING FORZADO Y EXCITADO PARAMETRICAMENTE
C                call duffing(NN,v,f,b,d,z1,e,pp,qq,fi,W)
                 f(1) = v(2)
                 f(2) = v(1)-b*v(1)**3-d*v(2)+z1*DCOS(v(3))
     *                 -b*e*(v(1)**3)*DCOS(pp*v(3)/qq+fi)
                 f(3) = W

                 do i=0,2
                    f(i+4) = v(i+7)
                    f(i+7) = (1-3*b*(v(1)**2)*(1+e*DCOS(pp*v(3)/qq+fi)))
     *                       *v(i+4)
     *                       -d*v(i+7)+((pp/qq)*b*e*(v(1)**3)
     *                       *DSIN(pp*v(3)/qq+fi))*v(i+10)
     *                       -z1*DSIN(v(3))*v(i+10)
                    f(i+10) = 0.
                 end do
              else
                 icount=icount+1
                 if (mod(icount,iorb).eq.0) then
C                   call f_clambda(N,NN,x,v,znorm,gsc,cum,clambda)
                    znorm(1)=0.0
                    do j=1,n
                       znorm(1)=znorm(1)+v(n*j+1)**2
                    end do
                    znorm(1)=sqrt(znorm(1))

                    do j=1,n
                       v(n*j+1)=v(n*j+1)/znorm(1)
                    end do

                    do j=2,n
                       do k=1,(j-1)
                          gsc(k)=0.
                          do l=1,n
                             gsc(k)=gsc(k)+v(n*l+j)*v(n*l+k)
                          end do
                       end do

                       do k=1,n
                         do l=1,(j-1)
                            v(n*k+j)=v(n*k+j)-gsc(l)*v(n*k+l)
                          end do
                       end do

                       znorm(j)=0.
                       do k=1,n
                          znorm(j)=znorm(j)+v(n*k+j)**2
                       end do
                       znorm(j)=sqrt(znorm(j))

                       do k=1,n
                          v(n*k+j)=v(n*k+j)/znorm(j)
                       end do
                    end do

                    do k=1,n
                       cum(k)=cum(k)+log(znorm(k))
                       clambda(k)=cum(k)/x
                    end do

                    if (icount.eq.nstep) then
C                     Escribe por pantalla
                      write(*,204)e,fi,clambda(1),clambda(2)
C                     Escribe a fichero lo mismo que antes
                      write(5,204)e,fi,clambda(1),clambda(2)
                      icount=0
                      loop=1
C                   end de icount.eq.nstep
                    endif
C                 end MOD(icount,iorb)Y
                  endif
C             end k.eq.1
              endif
C           END of WHILE
            END DO 
C        end jjj
         end do
C     end iii
      end do
      call cpu_time (time_end)
      PRINT *, 'Time of operation was ',time_end-time_begin,' seconds'
      close(5)
  204 format(2f14.5,2e16.5)
C     end program
      end

      function runge (n, y, f, x, h )
      implicit real*8(a-h,o-z)
      integer runge
      dimension phi(50),savey(50), y(n), f(n)
      data m/0/
      m=m+1
      go to (1,2,3,4,5),m

C      select case m  Se puede utilizar en F90/F95 no anterior y la estructura seria
C            case (1) 
Cc                paso 1
C            case (2)
Cc                paso 2
C            case (3)
Cc                paso 3
C            case (4)
Cc                paso 4
C            case (5)
Cc                paso 5
C      end select
c     paso 1
  1   runge=1
      return
c     paso 2
  2   do j=1,n
         savey(j)=y(j)
         phi(j)=f(j)
         y(j)=savey(j)+0.5*h*f(j)
      end do
      x=x+0.5*h
      runge=1
      return
c     paso 3
  3   do j=1,n
         phi(j)=phi(j)+2.0*f(j)
         y(j)=savey(j)+0.5*h*f(j)
      end do
      runge=1
      return
c     paso 4
  4   do j=1,n
         phi(j)=phi(j)+2.0*f(j)
         y(j)=savey(j)+h*f(j)
      end do
      x=x+0.5*h
      runge=1
      return
c     paso 5
  5   do j=1,n
         y(j)=savey(j)+(phi(j)+f(j))*h/6.0
      end do
      m=0.
      runge=0.
      return
      end

