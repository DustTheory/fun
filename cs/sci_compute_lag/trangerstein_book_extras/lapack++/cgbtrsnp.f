c modified from cgbtrs to avoid pivoting
      subroutine cgbtrsnp(trans,n,kl,ku,nrhs,ab,ldab,b,ldb,info)
      character trans
      integer info,kl,ku,ldab,ldb,n,nrhs
      complex ab(ldab,*),b(ldb,*)
      complex one
      parameter (one = (1.0,0.0))
      logical lnoti,notran
      integer i,j,kd,lm
      logical lsame
      external lsame
      external xerbla,cgemv,cgeru,clacgv,ctbsv
      intrinsic max,min

      info=0
      notran=lsame(trans,'N')
      if (.not.notran .and. .not.lsame(trans,'t') .and.
     &.not.lsame(trans,'C')) then
        info=-1
      else if (n.lt.0) then
        info=-2
      else if (kl.lt.0) then
        info=-3
      else if (ku.lt.0) then
        info=-4
      else if (nrhs.lt.0) then
        info=-5
      else if (ldab.lt.(kl+ku+1)) then
        info=-7
      else if (ldb.lt.max(1,n)) then
        info=-9
      endif
      if (info.ne.0) then
        call xerbla('cgbtrsnp',-info)
        return
      endif
      if (n.eq.0 .or. nrhs.eq.0) return
      kd=ku+1
      lnoti=kl.gt.0
      if (notran) then
        if (lnoti) then
          do j=1,n-1
            lm=min(kl,n-j)
            call cgeru(lm,nrhs,-one,ab(kd+1,j),1,b(j,1),ldb,
     &        b(j+1,1),ldb)
          enddo
        endif
        do i=1,nrhs
          call ctbsv('Upper','No transpose','Non-unit',n,ku,
     &      ab,ldab,b(1,i),1)
        enddo
      else if (lsame(trans,'T')) then
        do i=1,nrhs
          call ctbsv('Upper','Transpose','Non-unit',n,ku,ab,
     &      ldab,b(1,i),1)
        enddo
        if (lnoti) then
          do j=n-1,1,-1
            lm=min(kl,n-j)
            call cgemv('Transpose',lm,nrhs,-one,b(j+1,1),
     &        ldb,ab(kd+1,j),1,one,b(j,1),ldb)
          enddo
        endif
      else
        do i=1,nrhs
          call ctbsv('Upper','Conjugate transpose','Non-unit',n,
     &      ku,ab,ldab,b(1,i),1)
        enddo
        if (lnoti) then
          do j=n-1,1,-1
            lm=min(kl,n-j)
            call clacgv(nrhs,b(j,1),ldb)
            call cgemv('Conjugate transpose',lm,nrhs,-one,b(j+1,1),ldb,
     &        ab(kd+1,j),1,one,b(j,1),ldb)
            call clacgv(nrhs,b(j,1),ldb)
          enddo
        endif
      endif
      return
      end
