c from Dennis and Schnabel Algorithm A5.5.2
c replace A on and below diagonal with Cholesky factor L
c add mu to diagonal of A in order to keep entries of L < maxoffl 
      subroutine choldecomp(n,A,maxoffl,mu)
      integer k,n
      double precision A(n,n),maxoffl,mu
      
      double precision mu,minl,minl2,sqrtrnd
      integer j
   
      double precision roundoff
      common/const/roundoff

      sqrtrnd=sqrt(roundoff)
      minl=maxoffl * sqrt(sqrtrnd)
      if (maxoffl.le.0.d0) then
c       find max diagonal entry of A
        k=idamax(n,A(1,1),n+1)
        maxoffl=abs(A(k,k))
      endif
      minl2=maxoffl*sqrtrnd
   
      mu=0.d0
      do k=1,n
        i=idamax(n-k,A(k+1,k),1)+k
        minlkk=max(abs(A(i,k))/maxoffl,minl)
        if (A(k,k)>minlkk) then
          A(k,k)=sqrt(A(k,k))
        else
          minlkk=max(minlkk,minl2)
          mu=max(mu,minlkk**2 - A(k,k)) 
          A(k,k)=minlkk
        endif
        call dscal_simple(n-k,1.d0/A(k,k),A(k+1,k))
        do j=k+1,n 
          call daxpy_simple(n-j,A(j,k),A(j,k),A(j,j))
        enddo   
      enddo
   
      return
      end
