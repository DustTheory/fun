      subroutine fcn(n,x,fvec,fjac,ldfjac,iflag)
      implicit none
      integer(kind=4), intent(in) :: n,ldfjac,iflag
      real(kind=8), dimension(n), intent(in) :: x
      real(kind=8), dimension(n), intent(out) :: fvec
      real(kind=8), dimension(ldfjac,n), intent(out) :: fjac
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
      integer(kind=4) :: info
      real(kind=8), dimension(2) :: x,fvec
      real(kind=8), dimension(2,2) :: fjac
      real(kind=8), dimension(15) :: wa
      real(kind=8) :: tol
      real(kind=4) :: elapsed,reslt,start
      real(kind=4), dimension(2) :: tarray
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
