      subroutine tridiagonal_factor(n,A)
      int i,ip1,n
      real A(n,-1:1)
      do i=1,n-1
        ip1=i+1
        A(ip1,-1) = A(ip1,-1) / A(i,0)
        A(ip1,0) = A(ip1,0) - A(ip1,-1) * A(i,1)
      enddo
      return
      end
