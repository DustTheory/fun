module f90sub
contains
      subroutine sub(vector,matrix)
      real(kind=8), allocatable, dimension(:,:) :: matrix(:,:)
      real(kind=8), allocatable, dimension(:) :: vector
      integer i,j
        
      print *, "in sub"
      print *, "vector = ", &
     &  (vector(j),j=lbound(vector,1),ubound(vector,1))
      do i=lbound(matrix,1),ubound(matrix,1)
        print *, "matrix = ", (matrix(i,j), &
     &  j=lbound(matrix,2),ubound(matrix,2))
      enddo
      return
      end subroutine
end module f90sub
