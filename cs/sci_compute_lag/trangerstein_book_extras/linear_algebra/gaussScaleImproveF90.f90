      program main
      integer :: n,info
      character equed
      real(kind=8) :: rcond,ferr(1),berr(1)
      real(kind=8), allocatable :: A(:,:),LU(:,:),b(:),x(:)
      real(kind=8), allocatable :: r(:),c(:),work(:)
      integer, allocatable :: ipiv(:),iwork(:)

      n=1024
      allocate(A(n,n))
      allocate(LU(n,n))
      call random_number(A)
      LU=A

      allocate(b(n))
      call random_number(b)

      allocate(x(n))
      allocate(r(n))
      allocate(c(n))
      allocate(work(4*n))
      allocate(ipiv(n))
      allocate(iwork(n))
      call dgesvx('E','N',n,1,A,n,LU,n,ipiv,equed,r,c,b,n,x,n, &
     & rcond,ferr,berr,work,iwork,info)
      print *, "scale and improve: info,equed,rcond,ferr,berr = ", &
     & info,equed,rcond,ferr(1),berr(1)

      deallocate(iwork)
      deallocate(ipiv)
      deallocate(work)
      deallocate(c)
      deallocate(r)
      deallocate(x)
      deallocate(b)
      deallocate(LU)
      deallocate(A)

      stop
      end
