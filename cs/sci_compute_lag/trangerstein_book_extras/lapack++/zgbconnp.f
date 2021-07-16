c modified from zgbcon to avoid pivoting
      subroutine zgbconnp(norm,n,kl,ku,ab,ldab,anorm,rcond,
     &  work,rwork,info)
      character norm
      integer info,kl,ku,ldab,n
      double precision anorm,rcond
      double precision rwork(*)
      complex*16 ab(ldab,*),work(*)
      double precision one,zero
      parameter (one=1.0d+0,zero=0.0d+0)
      logical lnoti,onenrm
      character normin
      integer ix,j,kase,kase1,kd,lm
      double precision ainvnm,scale,smlnum
      complex*16 t,zdum
      integer isave(3)
      logical lsame
      integer izamax
      double precision dlamch
      complex*16 zdotc
      external lsame,izamax,dlamch,zdotc
      external xerbla,zaxpy,zdrscl,zlacn2,zlatbs
      intrinsic abs,dble,dimag,min
      double precision cabs1
      cabs1(zdum)=abs(dble(zdum))+abs(dimag(zdum))

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
        call xerbla('zgbconnp',-info)
        return
      endif
      rcond=zero
      if (n.eq.0) then
        rcond=one
        return
      else if (anorm.eq.zero) then
        return
      endif

      smlnum=dlamch('Safe minimum')
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
        call zlacn2(n,work(n+1),work,ainvnm,kase,isave)
        if (kase.ne.0) then
          if (kase.eq.kase1) then
            if (lnoti) then
              do j=1,n-1
                lm=min(kl,n-j)
                t=work(j)
                call zaxpy(lm,-t,ab(kd+1,j),1,work(j+1),1)
              enddo
            endif
            call zlatbs('Upper','No transpose','Non-unit',normin,n,
     &        ku,ab,ldab,work,scale,rwork,info)
          else
            call zlatbs('Upper','Conjugate transpose','Non-unit',
     &        normin,n,ku,ab,ldab,work,scale,rwork,info)
            if (lnoti) then
              do j=n-1,1,-1
                lm=min(kl,n-j)
                work(j)=work(j)-zdotc(lm,ab(kd+1,j),1,work(j+1),1)
              enddo
            endif
          endif
          normin='Y'
          if (scale.ne.one) then
            ix=izamax(n,work,1)
            if (scale.lt.cabs1(work(ix))*smlnum .or. scale.eq.zero)
     &         goto 40
            call zdrscl(n,scale,work,1)
          endif
          goto 10
        endif
      if (ainvnm.ne.zero) rcond=(one/ainvnm)/anorm
   40 continue
      return
      end
