c modified from sgbcon to avoid pivoting
      subroutine sgbconnp(norm,n,kl,ku,ab,ldab,anorm,rcond,
     &  work,iwork,info)
      character norm
      integer info,kl,ku,ldab,n
      real anorm, rcond
      integer iwork(*)
      real ab(ldab,*),work(*)
      real one,zero
      parameter (one=1.0,zero=0.0)
      logical lnoti,onenrm
      character normin
      integer ix,j,kase,kase1,kd,lm
      real ainvnm,scale,smlnum,t
      integer isave(3)
      logical lsame
      integer isamax
      real sdot,slamch
      external lsame,isamax,sdot,slamch
      external saxpy,slacn2,slatbs,srscl,xerbla
      intrinsic abs,min

      info=0
      onenrm=norm.eq.'1' .or. lsame(norm,'O')
      if (.not.onenrm .and. .not.lsame( norm,'I')) then
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
        info=-7
      end if
      if (info.ne.0) then
        call xerbla('sgbconnp',-info)
        return
      endif
      rcond=zero
      if (n.eq.0) then
        rcond=one
        return
      else if (anorm.eq.zero) then
        return
      end if
      smlnum=slamch('Safe minimum')
      ainvnm=zero
      normin='N'
      if (onenrm) then
        kase1=1
      else
        kase1=2
      end if
      kd=ku+1
      lnoti=kl.gt.0
      kase=0
   10 continue
      call slacn2(n,work(n+1),work,iwork,ainvnm,kase,isave)
      if (kase.ne.0) then
        if (kase.eq.kase1) then
          if (lnoti) then
            do j=1,n-1
              lm=min(kl,n-j)
              t=work(j)
              call saxpy(lm,-t,ab(kd+1,j),1,work(j+1),1)
            enddo
          endif ! lnoti
          call slatbs('Upper','No transpose','Non-unit',normin,n,ku,
     &      ab,ldab,work,scale,work(2*n+1),info)
        else ! kase.ne.kase1
          call slatbs('Upper','Transpose','non-unit',normin,n,ku,
     &      ab,ldab,work,scale,work(2*n+1),info)
          if (lnoti) then
            do j=n-1,1,-1
              lm=min(kl,n-j)
              work(j)=work(j)-sdot(lm,ab(kd+1,j),1,work(j+1),1)
            enddo
          endif ! lnoti
        endif
        normin='Y'
        if (scale.ne.one) then
          ix=isamax(n,work,1)
          if (scale.lt.abs(work(ix))*smlnum .or. scale.eq.zero) goto 40
          call srscl(n,scale,work,1)
        endif
        goto 10
      endif
      if (ainvnm.ne.zero) rcond=(one/ainvnm)/anorm
   40 continue
      return
      end
