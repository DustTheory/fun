      subroutine stride1(m,n,A)
      integer m,n
      double precision A(m,n)
      integer i,j

      do j=1,n
        do i=1,m
          A(i,j)=1.d0
        enddo
      enddo

      return
      end

      subroutine stridem(m,n,A)
      integer m,n
      double precision A(m,n)
      integer i,j

      do i=1,m
        do j=1,n
          A(i,j)=1.d0
        enddo
      enddo

      return
      end

      subroutine strider(m,n,colperm,rowperm,A)
      integer m,n
      integer colperm(n),rowperm(m)
      double precision A(0:m-1,0:n-1)
      integer i,j

      do i=1,m
        do j=1,n
          A(rowperm(i),colperm(j))=1.d0
        enddo
      enddo

      return
      end
