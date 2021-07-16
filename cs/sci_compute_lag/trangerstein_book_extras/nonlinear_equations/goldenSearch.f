      program main
      integer i,nsteps
      real a,b,c,d,f,fa,fb,fc,fd,fx,tau,x

      f(x)=x*x-4.*sin(x)

      tau=2./(1.+sqrt(5.))

      print *, "enter nsteps"
      read *, nsteps
      print *, "enter a"
      read *, a
      fa=f(a)
      print *, "enter b"
      read *, b
      fb=f(b)
      c=a+tau*(b-a)
      fc=f(c)
      if (fc.gt.max(fa,fb)) then
        print *, "function not unimodal on (",a,",",b,")"
        print *, "a,f(a) = ",a,fa
        print *, "c,f(c) = ",c,fc
        print *, "b,f(b) = ",b,fb
        call abort()
      endif
      d=b-tau*(b-a)
      fd=f(d)
      if (fd.gt.max(fa,fb)) then
        print *, "function not unimodal on (",a,",",b,")"
        print *, "a,f(a) = ",a,fa
        print *, "d,f(d) = ",d,fd
        print *, "b,f(b) = ",b,fb
        call abort()
      endif

      print *, -3,a,fa
      print *, -2,b,fb
      print *, -1,c,fc
      print *, 0,d,fd
      do i=1,nsteps
        if (fd.gt.fc) then
          a=d
          fa=fd
          d=c
          fd=fc
          c=a+tau*(b-a)
          fc=f(c)
          print *, i,c,fc
          if (fc.gt.max(fa,fb)) then
            print *, "function not unimodal on (",a,",",b,")"
            print *, "a,f(a) = ",a,fa
            print *, "c,f(c) = ",c,fc
            print *, "b,f(b) = ",b,fb
            call abort()
          endif
        else
          b=c
          fb=fc
          c=d
          fc=fd
          d=b-tau*(b-a)
          fd=f(d)
          print *, i,d,fd
          if (fd.gt.max(fa,fb)) then
            print *, "function not unimodal on (",a,",",b,")"
            print *, "a,f(a) = ",a,fa
            print *, "d,f(d) = ",d,fd
            print *, "b,f(b) = ",b,fb
            call abort()
          endif
        endif
      enddo

      stop
      end
