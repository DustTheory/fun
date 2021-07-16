      program main
      real(kind=8) eps,lambda,rho
      integer pow,i

      print *, "find smallest number eps such that 1 + eps > 1:"
      eps=epsilon(1.d0)
      do i=1,2
        lambda=1+eps
        rho=(lambda-1.d0)/eps
        pow=nint(log(eps)/log(2.d0))
        print '("eps = 2^{",i3,"} = ",z16," , ((1.+eps)-1.)/eps = ", &
     &    f2.0)',pow,eps,rho
        eps=eps*0.5d0
      enddo
      print *

      print *, "find largest number eps such that 1 + eps = 1:"
      eps=epsilon(1.d0)
      do i=1,3
        lambda=1.+eps;
        rho=(lambda-1.)/eps;
        pow=nint(log(eps)/log(2.d0))
        print '("eps = 2^{",i3,"} = ",z16," , ((1.+eps)-1.)/eps = ", &
     &    f2.0)',pow,eps,rho
        eps=eps*0.5d0
      enddo
      print *

      print *, &
     &  "find largest number eps such that ( 1 / eps - 1 ) * eps = 1:"
      eps=epsilon(1.d0)
      do i=1,3
        lambda=1./eps;
        rho=(lambda-1.)*eps;
        pow=nint(log(eps)/log(2.d0))
        print '("eps = 2^{",i3,"} = ",z16," , (1./eps-1.)*eps = ", &
     &    f18.16)',pow,eps,rho
        eps=eps*0.5d0
      enddo
      print *

      print *, "find largest number eps such that", &
     &  " ( 1 - 1 / eps ) + 1 / eps = 0:"
      eps=epsilon(1.d0)
      do i=1,3
        lambda=1./eps;
        rho=(1.-lambda)+lambda;
        pow=nint(log(eps)/log(2.d0))
        print '("eps = 2^{",i3,"} = ",z16," , (1.-1./eps)+1./eps = ", &
     &    f2.0)',pow,eps,rho
        eps=eps*0.5d0
      enddo

      stop
      end
