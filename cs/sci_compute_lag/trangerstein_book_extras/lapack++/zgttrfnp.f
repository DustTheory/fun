c modified from zgttrf to avoid pivoting
      subroutine zgttrfnp(n,dl,d,du,info)
      integer info,n
      complex*16 d(*),dl(*),du(*)
      double precision zero
      parameter (zero=0.0d+0)
      integer i
      complex*16 fact,zdum
      external xerbla
      intrinsic abs,dble,dimag
      double precision cabs1
      cabs1(zdum)=abs(dble(zdum))+abs(dimag(zdum))

      info=0
      if (n.lt.0) then
        info=-1
        call xerbla('zgttrfnp',-info)
        return
      endif
      if ( n.eq.0 ) return
      do i=1,n-2
        if (cabs1(d(i)).ne.zero) then
          fact=dl(i)/d(i)
          dl(i)=fact
          d(i+1)=d(i+1)-fact*du(i)
        endif
      enddo
      if (n.gt.1) then
        i=n-1
        if (cabs1(d(i)).ne.zero) then
          fact=dl(i)/d(i)
          dl(i)=fact
          d(i+1)=d(i+1)-fact*du(i)
        endif
      endif
      do i=1, n
        if (cabs1(d(i)).eq.zero) then
          info=i
          goto 50
        endif
      enddo
   50 continue
      return
      end
