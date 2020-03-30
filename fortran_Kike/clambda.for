      subroutine f_clambda(n,nn,x,v,znorm,gsc,cum,clambda)

      integer n,nn
      integer j,k,l
      real*8 x
      real*8 znorm(*),gsc(*),cum(*),clambda(*)
      real*8 v(*)

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
	 cum(k)=cum(k)+log(znorm(k))
	 clambda(k)=cum(k)/x
  90  continue
      end
