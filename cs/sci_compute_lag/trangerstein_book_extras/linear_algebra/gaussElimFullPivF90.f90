      program main
      integer :: i,j,n,info
      real(kind=8) :: s
      real(kind=8), allocatable :: A(:,:),LU(:,:),b(:),x(:)
      integer, allocatable :: ipiv(:),jpiv(:)

      n=10
      allocate(A(n,n))
      allocate(LU(n,n))
      call random_number(A)
      LU=A

      allocate(ipiv(n))
      allocate(jpiv(n))
      call dgetc2(n,LU,n,ipiv,jpiv,info)
      print *, "info = ",info
      print *, "ipiv = ",(ipiv(i),i=1,n-1)
      print *, "jpiv = ",(jpiv(j),j=1,n-1)

      allocate(b(n))
      allocate(x(n));
      call random_number(b)
      x=b
      call dgesc2(n,LU,n,x,ipiv,jpiv,s)
      print *, "scale = ",s
      print *, "x = ",(x(j),j=1,n)

      b=b-matmul(A,x)
      print *, "error = ",(b(i),i=1,n)

      deallocate(x)
      deallocate(b)
      deallocate(ipiv)
      deallocate(jpiv)
      deallocate(LU)
      deallocate(A)

      stop
      end
