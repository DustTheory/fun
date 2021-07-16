c     modified from dgtsv to avoid row pivoting
      subroutine dgtsvnp(n,nrhs,dl,d,du,b,ldb,info)
      integer info,ldb,n,nrhs
      double precision b(ldb,*),d(*),dl(*),du(*)

      double precision zero
      parameter (zero=0.0d+0)
      integer i,j
      double precision fact
      external xerbla
 
      info = 0
      if (n.lt.0) then
        info=-1
      else if (nrhs.lt.0) then
        info=-2
      else if (ldb.lt.max(1,n)) then
        info=-7
      endif
      if (info.ne.0) then
        call xerbla('dgtsvnp',-info)
        return
      endif

      if (n.eq.0) return

      if (nrhs.eq.1) then
        do i=1,n-2
          if (d(i).ne.zero) then
            fact=dl(i)/d(i)
            d(i+1)=d(i+1)-fact*du(i)
            b(i+1,1)=b(i+1,1)-fact*b(i,1)
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
            b(i+1,1)=b(i+1,1)-fact*b(i,1)
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
              b(i+1,j)=b(i+1,j)-fact*b(i,j)
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
               b(i+1,j)=b(i+1,j)-fact*b(i,j)
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
          b(n,j)=b(n,j)/d(n)
          if (n.gt.1) b(n-1,j)=(b(n-1,j)-du(n-1)*b(n,j))/d(n-1)
          do i=n-2,1,-1
            b(i,j)=(b(i,j)-du(i)*b(i+1,j)-dl(i)*b(i+2,j))/d(i)
          enddo
          if (j.lt.nrhs) then
            j=j+1
            goto 70
          endif
      else
        do j=1,nrhs
          b(n,j)=b(n,j)/d(n)
          if (n.gt.1) b(n-1,j)=(b(n-1,j)-du(n-1)*b(n,j))/d(n-1)
          do i=n-2,1,-1
            b(i,j)=(b(i,j)-du(i)*b(i+1,j)-dl(i)*b(i+2,j))/d(i)
          enddo
        enddo
      endif
      return
      end
