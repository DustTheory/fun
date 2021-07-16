      program main
      integer i,nsteps
      double precision f,fx,fpx,s,scal,sold,x,xn,z,zeta
      double precision scalf,scalx

      f(x)=x*x-4.d0*sin(x)
      fp(x)=2.d0*x-4.d0*cos(x)
      scal(x)=max(x*x,4.d0*abs(sin(x)))
c     z=1.9337537628270212d0
      z=0.d0

c     f(x)=1.d16*x*x-4.d0*sin(1.d8*x)
c     fp(x)=2.d16*x-4.d8*cos(1.d8*x)
c     scal(x)=max(1.d16*x*x,4.d0*abs(sin(1.d8*x)))
c     z=1.9337537628270212d-8
c     z=0.d0

c     f(x)=atan(x)
c     fp(x)=1.d0/(1.d0+x*x)
c     scal(x)=abs(atan(x))
c     z=0.d0

c     f(x)=x**2 - 2.d0
c     fp(x)=2.d0 * x
c     scal(x)=2.d0
c     z=sqrt(2.d0)

      print *, "enter nsteps"
      read *, nsteps
      print *, "enter x"
      read *, x

      sold=1.d0
      scalx=abs(x)
      scalf=abs(f(x))
      do i=1,nsteps
        fx=f(x)
        fpx=fp(x)
        xn=x-fx/fpx
        zeta=max(abs(x),abs(xn))
        s=xn-x
        print *,i,abs(x-z),zeta,scal(x),abs(s)/zeta,abs(fx)/scal(x),
     &    abs(s)/abs(sold)
c       print *,i,abs(x-z),scalx,scalf,abs(s)/scalx,abs(fx)/scalf,
c    &    abs(s)/abs(sold)
        if (xn.eq.x) stop
        sold=s
        x=xn
      enddo

      stop
      end
