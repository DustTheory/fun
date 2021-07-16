      PrOgRaM callingLapackFromFortran ! Fortran is insensitive to case
     ! by default, Fortran assumes that
     !   variable names beginning with (a-h) or (o-z) are real, and
     !   variable names beginning with (i-n) are integers
      implicit none
     ! now Fortran does not assume variable types, so
     !   we will have to declare all variables
     ! as a result, the compiler will help us to catch typos after
     !   declarations

     ! Fortran 90 style declarations :
      integer :: error, i, j, k ! Fortran 90 style declaration
      double precision ::  zero ! 64-bit floating point
      real, dimension(2) :: tarray ! 32-bit float array
      double precision, allocatable :: A(:,:) ! can allocate memory later
      double precision, allocatable :: y(:)
     !other data types:
     !  real : s,f ! 32-bit floating point number
     !  complex : c ! 32-bit complex floating point number
     !  double complex : z ! 64-bit complex floating point number
     !other array declaration styles:
     !  double precision :: x(512),y(512),A(512,512),B(512,512)
     !  double precision, dimension(512) :: x,y
     !  double precision, dimension(512,512) :: A,B

     ! Fortran 66 or 77 style declarations :
      integer m,n
      double precision xj,x(512)
      real elapsed,reslt,start
     !real s,f
     !complex c
     !double complex z

     !Fortran 90 external function declaration:
        double precision, external :: ddot 
     !Fortran 66 or 77 external function declaration:
     !  double precision ddot 
     !  external ddot 

     !We will also call daxpy and dgemv below.
     !  We do not need to declare them as externals,
     !  since they are not functions : do not return a function value

      parameter (zero=0.d0) ! 64-bit; "zero=0." would be 32-bit

      m=512
      n=512
      allocate(A(m,n),stat=error) ! dynamic memory allocation
      if (error.ne.0) print *, "unable to allocate memory for A"
     !dynamic memory allocation is recommended for large arrays

      allocate(y(m))
     !  if ",stat=..." is omitted, program aborts if allocation fails

     !fill array A with random numbers
      call random_number(A) ! F90 passes A's dimensions secretly
      call random_number(x)
     !fill vector x with zeros
      y(:) = zero ! Matlab-like loop notation
     !A(1:512,:) = 1.d0 ! alternative array assignment form

     !the following code contains 5 different ways to compute A*x
     !  and determine the user time required

     !dtime will return
     !  tarray(1) = User time in seconds
     !  tarray(2) = System time in seconds
     !  reslt = Run time since start in seconds
     !dtime cannot measure time much less than 0.001 seconds
     !if you get zero time, put the computation in a loop 
     !  do count=1,4096
     !    whatever you were going to do
     !  enddo
     !and divide the elapsed time by 4096

      call dtime(tarray,reslt) 
      start = tarray(1) ! start the user clock
      do k=1,4096 ! power of 2
        do j=1,n
          xj=x(j)
          do i=1,m
            y(i)=y(i)+A(i,j)*x(j)
          enddo
        enddo
      enddo
      call dtime(tarray,reslt) 
      elapsed = (tarray(1) - start) / 4096.d0
      print *, "time for double loop j over i = ", elapsed

      call dtime(tarray,reslt) 
      start = tarray(1)
      do k=1,4096
        do j=1,n
          call daxpy(m,x(j),A(1,j),1,y,1) ! slice is j'th column of A 
     !    note that we did not need to copy
     !      the j'th column of A to a vector
        enddo
      enddo
      call dtime(tarray,reslt) 
      elapsed = ( tarray(1) - start ) / 4096.d0
      print *, "time for loop j over daxpy = ", elapsed

      call dtime(tarray,reslt) 
      start = tarray(1)
      do k=1,4096
        do i=1,m
          y(i)=ddot(n,A(i,1),m,x,1) ! slice is i'th row of A 
     !    note that we did not need to copy
     !      the i'th row of A to a vector
        enddo
      enddo
      call dtime(tarray,reslt) 
      elapsed = ( tarray(1) - start ) / 4096.d0
      print *, "time for loop i over ddot = ", elapsed

      call dtime(tarray,reslt) 
      start = tarray(1)
      do k=1,4096
        call dgemv('N',m,n,1.d0,A,m,x,1,0.d0,y,1)
      enddo
      call dtime(tarray,reslt) 
      elapsed = ( tarray(1) - start ) / 4096.d0
      print *, "time for dgemv = ", elapsed

      call dtime(tarray,reslt) 
      start = tarray(1)
      do k=1,8192 ! power of 2
        y=matmul(A,x) ! F95 intrinsic function
      enddo
      call dtime(tarray,reslt) 
      elapsed = ( tarray(1) - start ) / 8192.d0
      print *, "time for dgemv = ", elapsed

      deallocate(A)
      deallocate(y)
     !deallocate(x) ! bad: x was not allocated dynamically
      stop
      end
