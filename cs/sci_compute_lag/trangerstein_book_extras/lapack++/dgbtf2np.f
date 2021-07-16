c modified from dgbtf2 to avoid pivoting
      subroutine dgbtf2np(m,n,kl,ku,ab,ldab,info)
      integer info,kl,ku,ldab,m,n
      double precision ab(ldab,*)
      double precision one,zero
      parameter (one=1.0d+0,zero=0.0d+0)
      integer j,ju,km
      integer idamax
      external idamax
      external dger,dscal,dswap,xerbla
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
        call xerbla('dgbtf2np',-info)
        return
      end if
      if (m.eq.0 .or. n.eq.0) return
      ju=1
      do j=1,min(m,n)
        km=min(kl,m-j)
        if (ab(ku+j,j).ne.zero) then
          ju=max(ju,min(j+ku+j-1,n))
          if (km.gt.0) then
            call dscal(km,one/ab(ku+1,j),ab(ku+2,j),1)
            if (ju.gt.j)
     &        call dger(km,ju-j,-one,ab(ku+2,j),1,ab(ku,j+1),ldab-1,
     &          ab(ku+1,j+1),ldab-1)
          endif
        else
          if (info.eq.0) info=j
        end if
      enddo
      return
      end
