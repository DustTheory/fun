c store reflectors times b in b
      subroutine sormqr_simple(m,n,A,tau,b)
      integer j,k,m,n
      real A(m,n),b(m),rho,tau(n)
      do k=1,min(m,n)
        rho=A(k,k)
        A(k,k)=1.
        alpha=tau(k)*sdot_simple(m-k+1,A(k,k),b(k))
        call saxpy_simple(m-k+1,-alpha,A(k,k),b(k))
        A(k,k)=rho
      enddo
      return
      end
