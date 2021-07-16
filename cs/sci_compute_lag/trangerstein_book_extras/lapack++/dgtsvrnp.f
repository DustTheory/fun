c modified from dgtsv to solve for rows of b, rather than columns
      subroutine dgtsvrnp(n,nrhs,dl,d,du,b,ldb,info)
      integer info,ldb,n,nrhs
      double precision b(ldb,*),d(*),dl(*),du(*)
      double precision zero
      parameter (zero=0.0d+0)
      integer i,j
      double precision fact
      intrinsic abs,max
      external xerbla

      info=0
      if (n.lt.0) then
        info=-1
      else if (nrhs.lt.0) then
        info=-2
      else if (ldb.lt.max(1,nrhs)) then
        info=-7
      endif
      if (info.ne.0) then
        call xerbla('dgtsvr',-info)
        return
      endif

      if (n.eq.0) return
      if (nrhs.eq.1) then
        do i=1,n-2
          if (d(i).ne.zero) then
            fact=dl(i)/d(i)
            d(i+1)=d(i+1)-fact*du(i)
            b(1,i+1)=b(1,i+1)- fact*b(1,i)
          else
            info=i
            return
          endif
          dl(i)=zero
        enddo
        if (n.gt.1) then
          i=n-1
          if (d(i).ne.zero) then
            fact=dl(i)/d(i)
            d(i+1)=d(i+1)-fact*du(i)
            b(1,i+1)=b(1,i+1)-fact*b(1,i)
          else
            info=i
            return
          endif
        endif
        if (d(n).eq.zero) then
          info=n
          return
        endif
      else
        do i=1,n-2
          if (d(i).ne.zero) then
            fact=dl(i)/d(i)
            d(i+1)=d(i+1)-fact*du(i)
            do j=1,nrhs
              b(j,i+1)=b(j,i+1)-fact*b(j,i)
            enddo
          else
            info=i
            return
          endif
          dl(i)=zero
        enddo
        if (n.gt.1) then
          i=n-1
          if (d(i).ne.zero) then
            fact=dl(i)/d(i)
            d(i+1)=d(i+1)-fact*du(i)
            do j=1,nrhs
              b(j,i+1)=b(j,i+1)-fact*b(j,i)
            enddo
          else
            info=i
            return
          endif
        endif
        if (d(n).eq.zero) then
          info=n
          return
        endif
      endif
      if (nrhs.le.2) then
        j=1
   70   continue
          b(j,n)=b(j,n)/d(n)
          if (n.gt.1) b(j,n-1)=(b(j,n-1)-du(n-1)*b(j,n))/d(n-1)
          do i=n-2,1,-1
            b(j,i)=(b(j,i)-du(i)*b(j,i+1)-dl(i)*b(j,i+2))/d(i)
          enddo
         if (j.lt.nrhs) then
           j=j+1
           goto 70
         endif
      else
         do j=1,nrhs
            b(j,n)=b(j,n)/d(n)
            if (n.gt.1) b(j,n-1)=(b(j,n-1)-du(n-1)*b(j,n))/d(n-1)
            do i=n-2,1,-1
              b(j,i)=(b(j,i)-du(i)*b(j,i+1)-dl(i)*b(j,i+2))/d(i)
            enddo
        enddo
      endif
      return
      end
