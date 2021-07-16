c modified from dgbcon to avoid pivoting
      subroutine dgbconnp(norm,n,kl,ku,ab,ldab,anorm,rcond,
     &  work,iwork,info)
      character norm
      integer info,kl,ku,ldab,n
      double precision anorm, rcond
      integer iwork(*)
      double precision ab(ldab,*),work(*)
      double precision one,zero
      parameter (one=1.0d+0,zero=0.0d+0)
      logical lnoti,onenrm
      character normin
      integer ix,j,kase,kase1,kd,lm
      double precision ainvnm,scale,smlnum,t
      integer isave(3)
      logical lsame
      integer idamax
      double precision ddot,dlamch
      external lsame,idamax,ddot,dlamch
      external daxpy,dlacn2,dlatbs,drscl,xerbla
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
        call xerbla('dgbconnp',-info)
        return
      endif
      rcond=zero
      if (n.eq.0) then
        rcond=one
        return
      else if (anorm.eq.zero) then
        return
      end if
      smlnum=dlamch('Safe minimum')
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
      call dlacn2(n,work(n+1),work,iwork,ainvnm,kase,isave)
      if (kase.ne.0) then
        if (kase.eq.kase1) then
          if (lnoti) then
            do j=1,n-1
              lm=min(kl,n-j)
              t=work(j)
              call daxpy(lm,-t,ab(kd+1,j),1,work(j+1),1)
            enddo
          endif ! lnoti
          call dlatbs('Upper','No transpose','Non-unit',normin,n,ku,
     &      ab,ldab,work,scale,work(2*n+1),info)
        else ! kase.ne.kase1
          call dlatbs('Upper','Transpose','non-unit',normin,n,ku,
     &      ab,ldab,work,scale,work(2*n+1),info)
          if (lnoti) then
            do j=n-1,1,-1
              lm=min(kl,n-j)
              work(j)=work(j)-ddot(lm,ab(kd+1,j),1,work(j+1),1)
            enddo
          endif ! lnoti
        endif
        normin='Y'
        if (scale.ne.one) then
          ix=idamax(n,work,1)
          if (scale.lt.abs(work(ix))*smlnum .or. scale.eq.zero) goto 40
          call drscl(n,scale,work,1)
        endif
        goto 10
      endif
      if (ainvnm.ne.zero) rcond=(one/ainvnm)/anorm
   40 continue
      return
      end
