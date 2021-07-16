      program main
      use f90sub
      real(kind=8), allocatable, dimension(:,:) :: matrix
      real(kind=8), allocatable, dimension(:) :: vector
      integer :: i,j,m,n

      m=3
      n=2
      allocate(matrix(m,n))
      allocate(vector(n))

      print *,"in main"
      do j=1,n
        do i=1,m
          matrix(i,j)=1.d0/dble(i+j-1)
        enddo
        vector(j)=1.d0
      enddo
      call sub(vector,matrix)
      print *,"back in main"

      stop
      end
