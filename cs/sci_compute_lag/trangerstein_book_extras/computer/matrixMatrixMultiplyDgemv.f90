      program matrixMatrixMultiply
      implicit none
      integer :: error,j,n ! Fortran 90 style declaration
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
      do j=1,n
        call dgemv('N',n,n,1.d0,A,n,B(1,j),1,0.d0,C(1,j),1)
      enddo
!     call dtime(tarray,reslt) 
!     print *, "time for loop j over dgemv = ", tarray(1) - start

      deallocate(C)
      deallocate(B)
      deallocate(A)
      stop
      end
