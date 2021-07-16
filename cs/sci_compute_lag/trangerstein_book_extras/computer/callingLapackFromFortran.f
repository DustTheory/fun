      PrOgRaM callingLapackFromFortran ! Fortran is insensitive to case
c     by default, Fortran assumes that
c       variable names beginning with (a-h) or (o-z) are real, and
c       variable names beginning with (i-n) are integers
      implicit none
c     now Fortran does not assume variable types, so
c       we will have to declare all variables
c     as a result, the compiler will help us to catch typos after
c       declarations

c     Fortran 66 or 77 style declarations :
      integer i,j,k,m,n
      double precision xj,x(512),y(512),zero ! 64-bit floating point
      real*8 A(512,512) ! 64-bit floating point
      real elapsed,reslt,tarray(2),start ! 32-bit floating point 
c     other data types:
c       complex c
c       double complex z

      double precision ddot 
      external ddot 

c     We will also call daxpy and dgemv below.
c       We do not need to declare them as externals,
c       since they are not functions : do not return a function value

      parameter (zero=0.d0) ! 64-bit; "zero=0." would be 32-bit

      m=512
      n=512

c     fill array A with random numbers
      do j=1,n
        do i=1,m
          A(i,j)=rand()
        enddo
        x(j)=rand()
      enddo
      do i=1,m
        y(i)=0.d0
      enddo

c     dtime will return
c       tarray(1) = User time in seconds
c       tarray(2) = System time in seconds
c       reslt = Run time since start in seconds
c     dtime cannot measure time much less than 0.001 seconds
c     if you get zero time, put the computation in a loop 
c       do count=1,4096
c         whatever you were going to do
c       enddo
c     and divide the elapsed time by 1000
      call dtime(tarray,reslt) 
      start = tarray(1) ! start the user clock
      do k=1,4096
        do j=1,n
          xj=x(j)
          do i=1,m
            y(i)=y(i)+A(i,j)*x(j)
          enddo
        enddo
      enddo
      call dtime(tarray,reslt) 
      elapsed = tarray(1) - start
      print *, "time for double loop j over i = ",elapsed / 4096.d0

      call dtime(tarray,reslt) 
      start = tarray(1)
      do k=1,4096
        do j=1,n
          call daxpy(m,x(j),A(1,j),1,y,1) ! slice is j'th column of A 
c         note that we did not need to copy
c           the j'th column of A to a vector
        enddo
      enddo
      call dtime(tarray,reslt) 
      elapsed = tarray(1) - start
      print *, "time for loop j over daxpy = ",elapsed / 4096.d0

      call dtime(tarray,reslt) 
      start = tarray(1)
      do k=1,4096
        do i=1,m
          y(i)=ddot(n,A(i,1),m,x,1) ! slice is i'th row of A 
c         note that we did not need to copy
c           the i'th row of A to a vector
        enddo
      enddo
      call dtime(tarray,reslt) 
      elapsed = tarray(1) - start
      print *, "time for loop i over ddot = ", elapsed / 4096.d0

      call dtime(tarray,reslt) 
      start = tarray(1)
      do k=1,4096
        call dgemv('N',m,n,1.d0,A,m,x,1,0.d0,y,1)
      enddo
      call dtime(tarray,reslt) 
      elapsed = tarray(1) - start
      print *, "time for dgemv = ", elapsed / 4096.d0

      stop
      end
