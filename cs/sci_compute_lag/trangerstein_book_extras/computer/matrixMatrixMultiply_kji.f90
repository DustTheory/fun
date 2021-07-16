      program matrixMatrixMultiply
      implicit none
      integer :: error,i,j,k,n ! Fortran 90 style declaration
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
      do k=1,n
        do j=1,n
          do i=1,n
            C(i,j)=C(i,j)+A(i,k)*B(k,j)
          enddo
        enddo
      enddo
!     call dtime(tarray,reslt) 
!     print *, "time for loop k over j over i = ", tarray(1) - start

      deallocate(C)
      deallocate(B)
      deallocate(A)
      stop
      end
