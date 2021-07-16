      subroutine fcn(n,x,fvec,fjac,ldfjac,iflag)
      implicit none
      integer n,ldfjac,iflag
      double precision x(n),fvec(n),fjac(ldfjac,n)

      if (iflag.eq.1) then
        fvec(1)=10.d0*(x(2)-x(1)**2)
        fvec(2)=1.d0-x(1)
      else
        fjac(1,1)=-20.d0*x(1)
        fjac(2,1)=-1.d0
        fjac(1,2)=10.d0
        fjac(2,2)=0.d0
      endif

      return
      end
      
      program callingMinpackFromFortran
      implicit none
      integer info
      double precision x(2),fvec(2),fjac(2,2),tol,wa(15)
      real elapsed,reslt,start,tarray(2)
      integer k
      external fcn

      tol=epsilon(1.d0)
      call dtime(tarray,reslt) 
      start=tarray(1)
      do k=1,4096
        x(1)=-1.d0
        x(2)=1.2d0
        call hybrj1(fcn,2,x,fvec,fjac,2,tol,info,wa,15)
      enddo
      call dtime(tarray,reslt) 
      elapsed=tarray(1)-start
      print *, "x = ",(x(k),k=1,2)
      print *, "average time for minpack hybrj1 = ",elapsed / 4096.d0

      stop
      end
