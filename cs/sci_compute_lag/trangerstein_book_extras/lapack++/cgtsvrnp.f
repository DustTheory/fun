c modified from cgtsvr to avoid pivoting
      subroutine cgtsvrnp(n,nrhs,dl,d,du,b,ldb,info)
      integer info,ldb,n,nrhs
      complex b(ldb,*),d(*),dl(*),du(*)
      complex zero
      parameter (zero=(0.0,0.0))
      integer j,k
      complex mult,zdum
      external xerbla
      real cabs1
      cabs1(zdum)=abs(real(zdum))+abs(aimag(zdum))

      info=0
      if (n.lt.0) then
        info=-1
      else if (nrhs.lt.0) then
        info=-2
      else if (ldb.lt.max(1,nrhs)) then
        info=-7
      end if
      if (info.ne.0) then
        call xerbla('cgtsvrnp',-info)
        return
      end if
      if (n.eq.0) return
      do k=1,n-1
        if (dl(k).eq.zero) then
          if (d(k).eq.zero) then
            info=k
            return
          end if
          mult=dl(k)/d(k)
          d(k+1)=d(k+1)-mult*du(k)
          do j=1,nrhs
            b(j,k+1)=b(j,k+1)-mult*b(j,k)
          enddo
          if (k.lt.(n-1)) dl(k)=zero
        end if
      enddo
      if (d(n).eq.zero) then
        info=n
        return
      end if
      do j=1,nrhs
        b(j,n)=b(j,n)/d(n)
        if (n.gt.1) b(j,n-1)=(b(j,n-1)-du(n-1)*b(j,n))/d(n-1)
        do k=n-2,1,-1
          b(j,k)=(b(j,k)-du(k)*b(j,k+1)-dl(k)*b(j,k+2))/d(k)
        enddo
      enddo
      return
      end
