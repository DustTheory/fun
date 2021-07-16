      program matrixMatrixMultiply
      implicit none
      integer :: error,i,j,k,n ! Fortran 90 style declaration
      real, dimension(2) :: tarray ! 32-bit float array
      real reslt,start
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

      call dtime(tarray,reslt) 
      start = tarray(1)
      do i=1,n
        do j=1,n
          do k=1,n
            C(i,j)=C(i,j)+A(i,k)*B(k,j)
          enddo
        enddo
      enddo
      call dtime(tarray,reslt) 
      print *, "time for loop i over j over k = ", tarray(1) - start

      call dtime(tarray,reslt) 
      start = tarray(1)
      do i=1,n
        do k=1,n
          do j=1,n
            C(i,j)=C(i,j)+A(i,k)*B(k,j)
          enddo
        enddo
      enddo
      call dtime(tarray,reslt) 
      print *, "time for loop i over k over j = ", tarray(1) - start

      call dtime(tarray,reslt) 
      start = tarray(1)
      do j=1,n
        do i=1,n
          do k=1,n
            C(i,j)=C(i,j)+A(i,k)*B(k,j)
          enddo
        enddo
      enddo
      call dtime(tarray,reslt) 
      print *, "time for loop j over i over k = ", tarray(1) - start

      call dtime(tarray,reslt) 
      start = tarray(1)
      do j=1,n
        do k=1,n
          do i=1,n
            C(i,j)=C(i,j)+A(i,k)*B(k,j)
          enddo
        enddo
      enddo
      call dtime(tarray,reslt) 
      print *, "time for loop j over k over i = ", tarray(1) - start

      call dtime(tarray,reslt) 
      start = tarray(1)
      do k=1,n
        do i=1,n
          do j=1,n
            C(i,j)=C(i,j)+A(i,k)*B(k,j)
          enddo
        enddo
      enddo
      call dtime(tarray,reslt) 
      print *, "time for loop k over i over j = ", tarray(1) - start

      call dtime(tarray,reslt) 
      start = tarray(1)
      do k=1,n
        do j=1,n
          do i=1,n
            C(i,j)=C(i,j)+A(i,k)*B(k,j)
          enddo
        enddo
      enddo
      call dtime(tarray,reslt) 
      print *, "time for loop k over j over i = ", tarray(1) - start

      call dtime(tarray,reslt) 
      start = tarray(1)
      do i=1,n
        do j=1,n
          C(i,j)=ddot(n,A(i,1),n,B(1,j),1)
        enddo
      enddo
      call dtime(tarray,reslt) 
      print *, "time for loop i over j over ddot = ", tarray(1) - start

      call dtime(tarray,reslt) 
      start = tarray(1)
      do j=1,n
        do i=1,n
          C(i,j)=ddot(n,A(i,1),n,B(1,j),1)
        enddo
      enddo
      call dtime(tarray,reslt) 
      print *, "time for loop j over i over ddot = ", tarray(1) - start

      call dtime(tarray,reslt) 
      start = tarray(1)
      do j=1,n
        do k=1,n
          call daxpy(n,B(k,j),A(1,k),1,C(1,j),1)
        enddo
      enddo
      call dtime(tarray,reslt) 
      print *, "time for loop j over k over daxpy = ", tarray(1) - start

      call dtime(tarray,reslt) 
      start = tarray(1)
      do k=1,n
        do j=1,n
          call daxpy(n,B(k,j),A(1,k),1,C(1,j),1)
        enddo
      enddo
      call dtime(tarray,reslt) 
      print *, "time for loop k over j over daxpy = ", tarray(1) - start

      call dtime(tarray,reslt) 
      start = tarray(1)
      do j=1,n
        call dgemv('N',n,n,1.d0,A,n,B(1,j),1,0.d0,C(1,j),1)
      enddo
      call dtime(tarray,reslt) 
      print *, "time for loop j over dgemv = ", tarray(1) - start

      call dtime(tarray,reslt) 
      start = tarray(1)
      call dgemm('N','N',n,n,n,1.d0,A,n,B,n,0.d0,C,n)
      call dtime(tarray,reslt) 
      print *, "time for dgemv = ", tarray(1) - start

      call dtime(tarray,reslt) 
      start = tarray(1)
      C=matmul(A,B) ! F95 intrinsic function
      call dtime(tarray,reslt) 
      print *, "time for matmul = ", tarray(1) - start

      deallocate(C)
      deallocate(B)
      deallocate(A)
      stop
      end
