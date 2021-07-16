c     modified from dgtcon to avoid pivoting
      subroutine dgtconnp(norm,n,dl,d,du,anorm,rcond,work,iwork,info)
      character norm
      integer info,n
      double precision anorm,rcond
      integer iwork(*)
      double precision d(*),dl(*),du(*),work(*)
      double precision one,zero
      parameter (one=1.0d+0,zero=0.0d+0)
      logical onenrm
      integer i,kase,kase1
      double precision ainvnm
      integer isave(3)
      logical lsame
      external lsame
      external dgttrs,dlacn2,xerbla

      info=0
      onenrm=norm.eq.'1' .or. lsame(norm,'O')
      if (.not.onenrm .and. .not.lsame(norm,'I')) then
        info=-1
      else if (n.lt.0) then
        info=-2
      else if (anorm.lt.zero) then
        info=-8
      endif
      if (info.ne.0) then
        call xerbla('dgtconnp',-info )
        return
      endif
      rcond=zero
      if (n.eq.0) then
        rcond=one
        return
      else if (anorm.eq.zero) then
        return
      end if
      do i=1,n
        if (d(i).eq.zero) return
      enddo
*
      ainvnm=zero
      if (onenrm) then
        kase1=1
      else
        kase1=2
      endif
      kase=0
   20 continue
        call dlacn2(n,work(n+1),work,iwork,ainvnm,kase,isave)
        if (kase.ne.0) then
        if (kase.eq.kase1) then
          call dgttrs('No transpose',n,1,dl,d,du,work,n,info)
        else
          call dgttrs('Transpose',n,1,dl,d,du,work,n,info)
        end if
        goto 20
      endif
      if (ainvnm.ne.zero) rcond=(one/ainvnm)/anorm
      return
      end
