c from Dennis and Schnabel Algorithm A5.5.1
c replace A on and below diagonal with Cholesky factor L
c add mu to diagonal of A in order to keep entries of L < maxoffl
      subroutine modelhess(n,A)
      integer n
      double precision A(n,n)

      integer k
      double precision maxoffl,mu,sqrtrnd

      double precision roundoff
      common/const/roundoff

      sqrtrnd=sqrt(roundoff)
      maxdiag=0.d0
      mindiag=0.d0
      do k=1,n
        absa=abs(A(k,k))
        maxdiag=max(maxdiag,absa)
        mindiag=min(mindiag,absa)
      enddo
      maxposdiag=max(0.d0,maxdiag)
      mu=0.d0
      if (mindiag.le.sqrtrnd*maxposdiag) then
        mu=2.d0*(maxposdiag-mindiag)*sqrtrnd-mindiag
        maxdiag=maxdiag+mu
      endif
      k=idamax(n*n,A(1,1),1)-1
      maxoff=abs(A(mod(k,n)+1,k/n+1))
      if (maxoff*(1.d0+2.d0*sqrtrnd).gt.maxdiag) then
        mu=mu+(maxoff-maxdiag)+2.d0*sqrtrnd*maxoff
        maxdiag=maxoff*(1.d0+2.d0*sqrtrnd)
      endif

      if (maxdiag.le.0.d0) then
        mu=1.d0
        maxdiag=1.d0
      endif

      if (mu.gt.0.d0) then
        do k=1,n
          A(k,k)=A(k,k)+mu
        enddo
      enddo
      maxoffl=sqrt(max(maxdiag,maxoff/dble(n)))
      call choldecomp(n,A,maxoffl,maxadd)

      if (maxadd.gt.0.d0) then
        maxev = A(1,1)
        minev = A(1,1)
        do k=1,n
          offrow=dasum(k-1,A(1,k),1)+dasum(n-k,A(k,k+1),n)
          maxev=max(maxev,A(k,k)+offrow)
          minev=min(minev,A(k,k)-offrow)
        enddo
        sdd=(maxev-minev)*sqrtrnd-minev
        sdd=max(sdd,0.d0)
        mu=min(maxadd,sdd)
        do k=1,n
          A(k,k)=A(k,k)+mu
        enddo
      endif
      call choldecomp(n,A,0.d0,maxadd)

      return
      end
