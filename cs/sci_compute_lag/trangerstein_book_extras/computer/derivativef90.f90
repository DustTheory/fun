      program main
      integer(kind=4) i,nlevels
      real(kind=8) diff,exact,f,fprime,h,x

! arithmetic statement functions:
!     f(x)=sqrt(x)
!     fprime(x)=0.5d0/sqrt(x)
      f(x)=sin(x)
      fprime(x)=cos(x)

      print *, "enter x"
      read *, x

      h=1.d0
      nlevels=64

      exact=fprime(x)
      do i=1,nlevels
        diff=(f(x+h)-f(x))/h
        write(11,*) -log10(h),log10(abs(diff-exact))
        h=h*0.5d0
      enddo

      stop
      end
