c     modified from sgttrs to avoid pivoting
      subroutine sgttrsnp(trans,n,nrhs,dl,d,du,b,ldb,info)
      character trans
      integer info,ldb,n,nrhs
      real b(ldb,*),d(*),dl(*),du(*)
      logical notran
      integer itrans,j,jb,nb
      integer ilaenv
      external ilaenv
      external sgtts2,xerbla
c     intrinsic max,min

      info=0
      notran=(trans.eq.'N' .or. trans.eq.'n')
      if (.not.notran .and. .not.(trans.eq.'T' .or. trans.eq.'t') .and.
     &  .not.(trans.eq.'C' .or. trans.eq.'c')) then
        info=-1
      else if (n.lt.0) then
        info=-2
      else if (nrhs.lt.0) then
        info=-3
      else if ( ldb.lt.max( n, 1 ) ) then
        info=-10
      end if
      if (info.ne.0) then
        call xerbla('sgttrsnp',-info)
        return
      end if
      if (n.eq.0 .or. nrhs.eq.0) return
      if (notran) then
        itrans=0
      else
        itrans=1
      end if
      if (nrhs.eq.1) then
        nb=1
      else
        nb=max(1,ilaenv(1,'sgttrs',trans,n,nrhs,-1,-1))
      end if

      if (nb.ge.nrhs) then
        call sgtts2(itrans,n,nrhs,dl,d,du,b,ldb)
      else
        do j=1,nrhs,nb
          jb=min(nrhs-j+1,nb)
          call sgtts2(itrans,n,jb,dl,d,du,b(1,j),ldb)
        enddo
      endif
      end
