      program main
      integer i,nsteps
      logical converged
      real a,b,delta,epsilon,f,fa,fb,fx,x

      f(x)=x*x-4.*sin(x)

      print *, "enter nsteps"
      read *, nsteps
      print *, "enter convergence tolerance on function"
      read *, epsilon
      print *, "enter convergence tolerance on solution"
      read *, delta
      print *, "enter lower bound on solution"
      read *, a
      fa=f(a)
      print *, -1,a,fa
      if (abs(fa).lt.epsilon) call exit()
      print *, "enter upper bound on solution"
      read *, b
      fb=f(b)
      print *, 0,b,fb
      if (abs(fb).lt.epsilon) call exit()
      if (abs(b-a).lt.delta) call exit()

      if (fa*fb.gt.0.) then
        print *, "zero not bracketed"
        call exit()
      endif

      i = 1
      converged=.false.
      do while (.not.converged .and. i.lt.nsteps)
        i=i+1
        x=0.5*(a+b)
        fx=f(x)
        converged=x.le.a .or. x.ge.b .or. abs(fx).lt.epsilon
        if (sign(1.,fx)*sign(1.,fa).gt.0.) then
          a=x
          fa=fx
        else
          b=x
          fb=fx
        endif
        converged=converged .or. abs(b-a).lt.delta
        print *, i,x,fx
      enddo

      stop
      end
