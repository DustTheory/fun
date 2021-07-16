c modified from sgbtf2 to avoid pivoting
      subroutine sgbtf2np(m,n,kl,ku,ab,ldab,info)
      integer info,kl,ku,ldab,m,n
      real ab(ldab,*)
      real one,zero
      parameter (one=1.0,zero=0.0)
      integer j,ju,km
      integer idamax
      external idamax
      external sger,sscal,sswap,xerbla
      intrinsic max,min

      info=0
      if (m.lt.0) then
        info=-1
      else if (n.lt.0) then
        info=-2
      else if (kl.lt.0) then
        info=-3
      else if (ku.lt.0) then
        info=-4
      else if (ldab.lt.kl+ku+1) then
        info=-6
      end if
      if (info.ne.0) then
        call xerbla('sgbtf2np',-info)
        return
      end if
      if (m.eq.0 .or. n.eq.0) return
      ju=1
      do j=1,min(m,n)
        km=min(kl,m-j)
        if (ab(ku+j,j).ne.zero) then
          ju=max(ju,min(j+ku+j-1,n))
          if (km.gt.0) then
            call sscal(km,one/ab(ku+1,j),ab(ku+2,j),1)
            if (ju.gt.j)
     &        call sger(km,ju-j,-one,ab(ku+2,j),1,ab(ku,j+1),ldab-1,
     &          ab(ku+1,j+1),ldab-1)
          endif
        else
          if (info.eq.0) info=j
        end if
      enddo
      return
      end
