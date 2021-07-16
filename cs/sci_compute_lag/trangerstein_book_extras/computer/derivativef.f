      program main
      integer i,nlevels
      double precision diff,exact,f,fprime,h,x

c arithmetic statement functions:
c     f(x)=sqrt(x)
c     fprime(x)=0.5d0/sqrt(x)
      f(x)=sin(x)
      fprime(x)=cos(x)

      print *, "enter x"
      read *, x

      h=1.d0
      nlevels=64

      exact=fprime(x)
c     print *, "exact = ",exact
      do i=1,nlevels
        diff=(f(x+h)-f(x))/h
        write(10,*) -log10(h),log10(abs(diff-exact))
        h=h*0.5d0
      enddo

      stop
      end
