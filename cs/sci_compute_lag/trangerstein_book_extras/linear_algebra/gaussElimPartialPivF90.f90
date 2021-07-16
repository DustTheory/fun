      program main
      integer :: i,j,n,info
      real(kind=8), allocatable :: A(:,:),LU(:,:),b(:),x(:)
      integer, allocatable :: ipiv(:)

      n=10
      allocate(A(n,n))
      allocate(LU(n,n))
      call random_number(A)
      LU=A

      allocate(ipiv(n))
      call dgetrf(n,n,LU,n,ipiv,info)
      print *, "info = ",info
      print *, "ipiv = ",(ipiv(i),i=1,n-1)

      allocate(b(n))
      allocate(x(n));
      call random_number(b)
      x=b
      call dgetrs('N',n,1,LU,n,ipiv,x,n,info)
      print *, "x = ",(x(j),j=1,n)

      b=b-matmul(A,x)
      print *, "error = ",(b(i),i=1,n)

      deallocate(x)
      deallocate(b)
      deallocate(ipiv)
      deallocate(LU)
      deallocate(A)

      stop
      end
