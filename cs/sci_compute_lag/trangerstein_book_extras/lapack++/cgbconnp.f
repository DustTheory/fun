c modified from cgbcon to avoid pivoting
      subroutine cgbconnp(norm,n,kl,ku,ab,ldab,anorm,rcond,
     &  work,rwork,info)
      character norm
      integer info,kl,ku,ldab,n
      real anorm,rcond
      real rwork(*)
      complex ab(ldab,*),work(*)
      real one,zero
      parameter (one=1.0,zero=0.0)
      logical lnoti,onenrm
      character normin
      integer ix,j,kase,kase1,kd,lm
      real ainvnm,scale,smlnum
      complex t,zdum
      integer isave(3)
      logical lsame
      integer izamax
      real slamch
      complex cdotc
      external lsame,izamax,slamch,cdotc
      external xerbla,caxpy,csrscl,clacn2,clatbs
      real cabs1
      cabs1(zdum)=abs(real(zdum))+abs(aimag(zdum))

      info=0
      onenrm=norm.eq.'1' .or. lsame(norm,'O')
      if (.not.onenrm .and. .not.lsame(norm,'I')) then
        info=-1
      else if (n.lt.0) then
        info=-2
      else if (kl.lt.0) then
        info=-3
      else if (ku.lt.0) then
        info=-4
      else if (ldab.lt.kl+ku+1) then
        info=-6
      else if (anorm.lt.zero) then
        info=-8
      endif
      if (info.ne.0) then
        call xerbla('cgbconnp',-info)
        return
      endif
      rcond=zero
      if (n.eq.0) then
        rcond=one
        return
      else if (anorm.eq.zero) then
        return
      endif

      smlnum=slamch('Safe minimum')
      ainvnm=zero
      normin='n'
      if (onenrm) then
        kase1=1
      else
        kase1=2
      endif
      kd=ku+1
      lnoti=kl.gt.0
      kase=0
   10 continue
        call clacn2(n,work(n+1),work,ainvnm,kase,isave)
        if (kase.ne.0) then
          if (kase.eq.kase1) then
            if (lnoti) then
              do j=1,n-1
                lm=min(kl,n-j)
                t=work(j)
                call caxpy(lm,-t,ab(kd+1,j),1,work(j+1),1)
              enddo
            endif
            call clatbs('Upper','No transpose','Non-unit',normin,n,
     &        ku,ab,ldab,work,scale,rwork,info)
          else
            call clatbs('Upper','Conjugate transpose','Non-unit',
     &        normin,n,ku,ab,ldab,work,scale,rwork,info)
            if (lnoti) then
              do j=n-1,1,-1
                lm=min(kl,n-j)
                work(j)=work(j)-cdotc(lm,ab(kd+1,j),1,work(j+1),1)
              enddo
            endif
          endif
          normin='Y'
          if (scale.ne.one) then
            ix=izamax(n,work,1)
            if (scale.lt.cabs1(work(ix))*smlnum .or. scale.eq.zero)
     &         goto 40
            call csrscl(n,scale,work,1)
          endif
          goto 10
        endif
      if (ainvnm.ne.zero) rcond=(one/ainvnm)/anorm
   40 continue
      return
      end
