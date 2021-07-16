      program main
      integer i,nsteps
      logical converged
      real a,b,c,f,fa,fb,fx,x,z

      f(x)=x+2.*(x-5.)
      g(x)=1.-2.*exp(-abs(f(x))/c**2)
      h(x)=(1./c)/((f(x)/c)+c**2/(f(x)/c))

c     print *, "enter nsteps"
c     read *, nsteps
c     print *, "enter a"
c     read *, a
c     print *, "enter b"
c     read *, b

      z=10./3.
      c=abs(f(10./3.))

      nsteps=30
      a=1.
      b=20./3.
      fa=f(a)
      fb=f(b)
      print *, -1,a,fa,g(a),h(a),1./fa
      print *, 0,b,fb,g(b),h(b),1./fb

      i = 1
      converged=.false.
      do while (.not.converged .and. i.lt.nsteps)
        i=i+1
        x=0.5*(a+b)
        fx=f(x)
        converged=x.le.a .or. x.ge.b
        if (fx*fa.gt.0.) then
          a=x
          fa=fx
        else
          b=x
          fb=fx
        endif
        print *, i,x,fx,g(x),h(x),1./fx
      enddo
c 3.*x-10 = 2**(-21)
      print *, 3.*x-10.
      print *, 3.*(x-5.)+5.

      stop
      end
