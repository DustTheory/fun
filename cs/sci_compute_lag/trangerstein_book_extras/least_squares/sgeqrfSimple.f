c store R on an above the diagonal of A
c store Householder vectors below the diagonal of A
      subroutine sgeqrf_simple(m,n,A,tau)
      integer j,k,m,n
      real A(m,n),prod,alpha,beta,tau(n),term
      do k=1,min(m,n)
        beta=-sign(snrm2(m-k+1,A(k,k),1),A(k,k))
        alpha=A(k,k)
        tau(k)=(beta-alpha)/beta
        call sscal(m-k,1./(alpha-beta),A(k+1,k),1)
        A(k,k)=one
        do j=k+1,n
          prod=sdot(m-k+1,A(k,k),A(k,j))
          call saxpy(m-k+1,-prod,A(k,k),A(k,j))
        enddo
        A(k,k)=beta
      enddo
      return
      end
