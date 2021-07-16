c     modified from dgttrf to avoid pivoting
      subroutine dgttrfnp(n,dl,d,du,info)
      integer info,n
      double precision d(*),dl(*),du(*)
      double precision zero
      parameter (zero=0.0d+0)
      integer i
      double precision fact
      intrinsic abs
      external xerbla

      info=0
      if (n.lt.0) then
        info=-1
        call xerbla('dgttrfnp',-info)
        return
      endif
      if (n.eq.0) return
      do i=1,n-2
        if (d(i).ne.zero) then
          fact=dl(i)/d(i)
          dl(i)=fact
          d(i+1)=d(i+1)-fact*du(i)
        endif
      enddo
      if (n.gt.1) then
        i=n-1
        if (d(i).ne.zero) then
          fact=dl(i)/d(i)
          dl(i)=fact
          d(i+1)=d(i+1)-fact*du(i)
        endif
      endif
      do i=1,n
        if (d(i).eq.zero) then
          info=i
          goto 50
        endif
      enddo
   50 continue
      return
      end
