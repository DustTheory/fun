c store R on an above the diagonal of A
c store Householder vectors below the diagonal of A
      subroutine sgeqrf_simple(m,n,A,tau,work)
      integer j,k,m,n
      real A(m,n),prod,alpha,beta,tau(n),term,work(n)
      do k=1,min(m,n)
        beta=-sign(snrm2(m-k+1,A(k,k),1),A(k,k))
        alpha=A(k,k)
        tau(k)=(beta-alpha)/beta
        call sscal(m-k,1./(alpha-beta),A(k+1,k),1)
        A(k,k)=one ! u_k now stored on and below diagonal of column k
        call sgemv('T',m-k,n-k-1,1.,A(k,k+1),m,A(k,k),1,0.,work,1)
        call sger(m-k,n-k-1,-tau(k),A(k,k),1,work,1,A(k,k+1),m)
        A(k,k)=beta ! since first entry of u_k = 1, can store beta
      enddo
      return
      end

c store reflectors times B in B
      subroutine sormqr_simple(m,n,k,A,tau,B)
      integer i,j,k,m,n
      real A(m,n),B(m,k),rho,tau(n)
      do i=1,min(m,n)
        rho=A(i,i)
        A(i,i)=1.
        call sgemv('T',m-i,k,1.,B(i,1),m,A(i,i),1,0.,work,1)
        call sger(m-i,k,-tau(i),A(i,i),1,work,1,B(i,1),m)
        A(i,i)=rho
      enddo
      return
      end
