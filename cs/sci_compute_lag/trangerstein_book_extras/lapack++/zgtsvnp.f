c modified from zgtsv to avoid pivoting
      subroutine zgtsvnp(n,nrhs,dl,d,du,b,ldb,info)
      integer info,ldb,n,nrhs
      complex*16 b(ldb,*),d(*),dl(*),du(*)
      complex*16 zero
      parameter (zero=(0.0d+0,0.0d+0))
      integer j,k
      complex*16 mult,zdum
      intrinsic abs,dble,dimag,max
      external xerbla
      double precision cabs1
      cabs1(zdum)=abs(dble(zdum))+abs(dimag(zdum))

      info=0
      if (n.lt.0) then
        info=-1
      else if (nrhs.lt.0) then
        info=-2
      else if (ldb.lt.max(1,n)) then
        info=-7
      endif
      if (info.ne.0) then
        call xerbla('zgtsvnp',-info)
        return
      endif
      if (n.eq.0) return
      do k=1, n - 1
        if (dl(k).eq.zero) then
          if (d(k).eq.zero) then
            info=k
            return
          endif
        endif
        mult=dl(k)/d(k)
        d(k+1)=d(k+1)-mult*du(k)
        do j=1,nrhs
          b(k+1,j)=b(k+1,j)-mult*b(k,j)
        enddo
        if (k.lt.(n-1)) dl(k)=zero
      enddo
      if (d(n).eq.zero) then
        info=n
        return
      endif
      do j=1, nrhs
        b(n,j)=b(n,j)/d(n)
        if (n.gt.1) b(n-1,j)=(b(n-1,j)-du(n-1)*b(n,j))/d(n-1)
        do k=n-2,1,-1
          b(k,j)=(b(k,j)-du(k)*b(k+1,j)-dl(k)*b(k+2,j))/d(k)
        enddo
      enddo
      return
      end
