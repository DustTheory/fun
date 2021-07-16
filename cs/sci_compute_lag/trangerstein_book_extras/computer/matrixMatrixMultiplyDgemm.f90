      program matrixMatrixMultiply
      implicit none
      integer :: error,n ! Fortran 90 style declaration
!     real, dimension(2) :: tarray ! 32-bit float array
!     real reslt,start
      double precision, allocatable :: A(:,:),B(:,:),C(:,:)
      double precision, external :: ddot 

      n=1024
      allocate(A(n,n),stat=error) ! dynamic memory allocation
      if (error.ne.0) print *, "unable to allocate memory for A"
     !dynamic memory allocation is recommended for large arrays

      allocate(B(n,n),stat=error) ! dynamic memory allocation
      if (error.ne.0) print *, "unable to allocate memory for B"

      allocate(C(n,n),stat=error) ! dynamic memory allocation
      if (error.ne.0) print *, "unable to allocate memory for C"

      call random_number(A)
      call random_number(B)
      C(:,:) = 0.d0 ! Matlab-like loop notation

!     call dtime(tarray,reslt) 
!     start = tarray(1)
      call dgemm('N','N',n,n,n,1.d0,A,n,B,n,0.d0,C,n)
!     call dtime(tarray,reslt) 
!     print *, "time for dgemv = ", tarray(1) - start

      deallocate(C)
      deallocate(B)
      deallocate(A)
      stop
      end
