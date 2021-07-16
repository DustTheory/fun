      subroutine tridiagonal_solve(n,A,b)
      int i,n
      real A(n,-1:1),b(n)
      do i=2,n
        b(i) = b(i) - A(i,-1) * b(i-1)
      enddo
      b(n) = b(n) / A(n,0)
      do i=n-1,1,-1
        b(i) = ( b(i) - A(i,1) * b(i+1) ) / A(i,0)
      enddo
      return
      end
