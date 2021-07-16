      program main
      integer m,n
      parameter (m=3,n=2)
      double precision matrix(m,n),vector(n)
      integer i,j

      print *,"in main"
      do j=1,n
        do i=1,m
          matrix(i,j)=1.d0/dble(i+j-1)
        enddo
        vector(j)=1.d0
      enddo
      call sub(m,n,vector,matrix)
      print *,"back in main"

      stop
      end
