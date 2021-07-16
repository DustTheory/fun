      subroutine zgtconnp(norm,n,dl,d,du,anorm,rcond,work,info)
      character norm
      integer info, n
      double precision anorm,rcond
      complex*16 d(*),dl(*),du(*),work(*)
      double precision one,zero
      parameter (one=1.0d+0,zero=0.0d+0)
      logical onenrm
      integer i,kase,kase1
      double precision ainvnm
      integer isave(3)
      logical lsame
      external lsame
      external xerbla,zgttrs,zlacn2
      intrinsic dcmplx

      info=0
      onenrm=norm.eq.'1' .or. lsame(norm,'O')
      if (.not.onenrm .and. .not.lsame(norm,'I')) then
        info=-1
      else if (n.lt.0) then
        info=-2
      else if (anorm.lt.zero) then
        info=-8
      endif
      if ( info.ne.0 ) then
        call xerbla( 'zgtconnp', -info )
        return
      endif
      rcond=zero
      if (n.eq.0) then
        rcond=one
        return
      else if (anorm.eq.zero) then
        return
      endif
      do i=1,n
        if (d(i).eq.dcmplx(zero)) return
      enddo
      ainvnm=zero
      if (onenrm ) then
        kase1=1
      else
        kase1=2
      endif
      kase=0
   20 continue
      call zlacn2(n,work(n+1),work,ainvnm,kase,isave)
      if (kase.ne.0) then
        if (kase.eq.kase1) then
          call zgttrsnp('No transpose',n,1,dl,d,du,work,n,info)
        else
          call zgttrsnp('Conjugate transpose',n,1,dl,d,du,work,n,info )
        endif
        goto 20
      endif
      if (ainvnm.ne.zero) rcond=(one/ainvnm)/anorm
      return
      end
