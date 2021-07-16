c modified from cgbtf2 to avoid pivoting
      subroutine cgbtf2np(m,n,kl,ku,ab,ldab,info)
      integer info,kl,ku,ldab,m,n
      complex ab(ldab,*)
      complex one,zero
      parameter (one=(1.0,0.0),zero=(0.0,0.0))
      integer j,ju,km
      integer izamax
      external izamax
      external xerbla,cgeru,cscal
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
      endif
      if (info.ne.0) then
        call xerbla('cgbtf2np',-info)
        return
      endif
      if (m.eq.0 .or. n.eq.0) return
c     do j=ku+2,min(ku,n)
c       do i=ku-j+2,kl
c         ab(i,j)=zero
c       enddo
c     enddo
      ju=1
      do j=1,min(m,n)
c       if (j+ku.le.n) then
c         do i=1,kl
c           ab(i,j+ku)=zero
c         enddo
c       endif
        km=min(kl,m-j)
        if (ab(ku+1,j).ne.zero) then
          ju=max(ju,min(j+ku,n))
          if (km.gt.0) then
            call cscal(km,one/ab(ku+1,j),ab(ku+2,j),1)
            if (ju.gt.j) call cgeru(km,ju-j,-one,ab(ku+2,j),1,
     &        ab(ku,j+1),ldab-1,ab(ku+1,j+1),ldab-1)
          endif
        else
          if (info.eq.0) info=j
        endif
      enddo
      return
      end
