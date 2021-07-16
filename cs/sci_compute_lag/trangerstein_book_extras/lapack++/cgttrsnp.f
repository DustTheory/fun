c modified from cgttrs to avoid pivoting
      subroutine cgttrsnp(trans,n,nrhs,dl,d,du,b,ldb,info)
      character trans
      integer info,ldb,n,nrhs
      complex b(ldb,*),d(*),dl(*),du(*)
      logical notran
      integer itrans,j,jb,nb
      integer ilaenv
      external ilaenv
      external xerbla,cgtts2
c     intrinsic max,min

      info=0
      notran=(trans.eq.'N' .or. trans.eq.'N')
      if (.not.notran .and. .not.(trans.eq.'T' .or. trans.eq.'t' )
     &.and. .not.( trans.eq.'C' .or. trans.eq.'c' ) ) then
        info=-1
      else if (n.lt.0) then
        info=-2
      else if (nrhs.lt.0) then
        info=-3
      else if (ldb.lt.max(n,1)) then
        info=-10
      endif
      if (info.ne.0) then
        call xerbla('cgttrsnp',-info)
        return
      endif
      if (n.eq.0 .or. nrhs.eq.0) return
      if (notran) then
        itrans=0
      else if (trans.eq.'T' .or. trans.eq.'t') then
        itrans=1
      else
        itrans=2
      endif
      if (nrhs.eq.1) then
        nb=1
      else
        nb=max(1, ilaenv(1,'cgttrs',trans,n,nrhs,-1,-1))
      endif
      if (nb.ge.nrhs) then
        call cgtts2np(itrans,n,nrhs,dl,d,du,b,ldb)
      else
        do j=1,nrhs,nb
          jb=min(nrhs-j+1,nb)
          call cgtts2np(itrans,n,jb,dl,d,du,b(1,j),ldb)
        enddo
      endif
      end
