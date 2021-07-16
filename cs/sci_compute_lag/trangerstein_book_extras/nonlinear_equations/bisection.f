      program main
      integer i,nsteps
      real a,b,f,fa,fb,fx,x

      f(x)=x*x-4.*sin(x)

      print *, "enter nsteps"
      read *, nsteps
      print *, "enter a"
      read *, a
      fa=f(a)
      print *, "enter b"
      read *, b
      fb=f(b)
      if (fa*fb.gt.0.) then
        print *, "zero not bracketed:"
        print *, "a,f(a) = ",a,fa
        print *, "b,f(b) = ",b,fb
        call exit()
      endif
      print *, -1,a,fa
      print *, 0,b,fb

      do i=1,nsteps
        x=0.5*(a+b)
        fx=f(x)
        if (fx*fa.gt.0.) then
          a=x
          fa=fx
        else
          b=x
          fb=fx
        endif
        print *, i,x,fx
      enddo

      stop
      end
