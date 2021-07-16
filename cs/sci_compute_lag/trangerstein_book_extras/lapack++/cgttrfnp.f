c modified from cgttrf to avoid pivoting
      subroutine cgttrfnp(n,dl,d,du,info)
      integer info,n
      complex d(*),dl(*),du(*)
      real zero
      parameter (zero=0.0)
      integer i
      complex fact,zdum
      external xerbla
      intrinsic sngl
      real cabs1
      cabs1(zdum)=abs(real(zdum))+abs(aimag(zdum))

      info=0
      if (n.lt.0) then
        info=-1
        call xerbla('cgttrfnp',-info)
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
