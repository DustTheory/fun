      program main
      integer :: i,j,m
      real, allocatable :: L(:,:),b(:),y(:),Ly(:)

      m=10
      allocate(L(m,m))
      do j=1,m
        L(j,j)=1.
        do i=j+1,m
          L(i,j)=rand()
        enddo
      enddo

      allocate(b(m))
      call random_number(b)

      allocate(y(m))
      y=b
      call strsv('L','N','U',m,L,m,y,1) !solve L*y = b, store soln in y
      print *, "y = ",(y(j),j=1,m)

      allocate(Ly(m))
      Ly=y
      call strmv('L','N','U',m,L,m,Ly,1) ! multiply Ly = L * y
!     DO NOT USE MATMUL TO COMPUTE L * y
!     MATMUL USES ENTRIES OF L ABOVE THE DIAGONAL
      b=b-Ly
      print *, "error = ",(b(i),i=1,m)

      deallocate(Ly)
      deallocate(y)
      deallocate(b)
      deallocate(L)
      stop
      end
