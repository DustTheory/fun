c modified from cgtrfs to work with rows of b and x
      subroutine cgtrfsr(trans,n,nrhs,dl,d,du,dlf,df,duf,du2,ipiv,b,ldb,
     &  x,ldx,ferr,berr,work,rwork,info)
      character trans
      integer info,ldb,ldx,n,nrhs
      integer ipiv(*)
      real berr(*),ferr(*),rwork(*)
      complex b(ldb,*),d(*),df(*),dl(*),dlf(*),du(*),du2(*),duf(*),
     &  work(*),x(ldx,*)
      integer itmax
      parameter (itmax=5)
      real zero,one
      parameter (zero=0.0,one=1.0)
      real two
      parameter (two=2.0)
      real three
      parameter (three=3.0)
      logical notran
      character ltrans,transn,transt
      integer count,i,j,kase,nz
      real eps,lstres,s,safe1,safe2,safmin
      complex zdum
      integer isave(3)
      external xerbla,zaxpy,zcopy,cgttrs,zlacn2,zlagtm
      logical lsame
      real dlamch
      external lsame,dlamch
      real cabs1
      cabs1(zdum)=abs(real(zdum))+abs(aimag(zdum))

      info=0
      notran=lsame(trans,'N')
      if (.not.notran .and. .not.lsame(trans,'T') .and.
     &.not.lsame(trans,'C')) then
        info=-1
      else if (n.lt.0) then
        info=-2
      else if (nrhs.lt.0) then
        info=-3
      else if (ldb.lt.max(1,nrhs)) then
        info=-13
      else if (ldx.lt.max(1,nrhs)) then
        info=-15
      endif
      if (info.ne.0) then
        call xerbla('cgtrfsr',-info)
        return
      endif
      if (n.eq.0 .or. nrhs.eq.0) then
        do j=1,nrhs
          ferr(j)=zero
          berr(j)=zero
        enddo
        return
      endif
      if (notran) then
        ltrans='T'
        transn='N'
        transt='C'
      else
        ltrans='N'
        transn='C'
        transt='N'
      endif
c         if trans==N: X A = B ==> A^T X^T = B^T ==> ltrans=T
c         if trans==T: X A^T = B ==> A X^T = B^T ==> ltrans=N
c         if trans==C: X A^H = B ==> A X^H = B^H ==> ltrans=N
      nz=4
      eps=dlamch('Epsilon')
      safmin=dlamch('Safe minimum')
      safe1=nz*safmin
      safe2=safe1/eps
      do j=1,nrhs
        count=1
        lstres=three
   20   continue
          call zcopy(n,b(j,1),ldb,work,1)
          if (lsame(trans,'C')) then
            do i=1,n
              work(i)=conjg(work(i))
              x(j,i)=conjg(x(j,i))
            enddo
          endif
c         call zlagtm(ltrans,n,1,-one,dl,d,du,x(j,1),ldx,one,work,n)
          if (notran) then
            call cgtmv(n,cmplx(-one),du,d,dl,x(j,1),ldx,cmplx(one),
     &        work,1)
          else
            call cgtmv(n,cmplx(-one),dl,d,du,x(j,1),ldx,cmplx(one),
     &        work,1)
          endif
          if (lsame(trans,'C')) then
            do i=1,n
              work(i)=conjg(work(i))
              x(j,i)=conjg(x(j,i))
            enddo
          endif
          if (notran) then
            if (n.eq.1) then
              rwork(1)=cabs1(b(j,1))+cabs1(d(1))*cabs1(x(j,1))
            else
              rwork(1)=cabs1(b(j,1))+cabs1(d(1))*cabs1(x(j,1))
     &          +cabs1(du(1))*cabs1(x(j,2))
              do i=2,n - 1
                rwork(i)=cabs1(b(j,i))+cabs1(dl(i-1))*cabs1(x(j,i-1))
     &            +cabs1(d(i))*cabs1(x(j,i))
     &            +cabs1(du(i))*cabs1(x(j,i+1))
              enddo
              rwork(n)=cabs1(b(j,n))+cabs1(dl(n-1))*cabs1(x(j,n-1))
     &          +cabs1(d(n))*cabs1(x(j,n))
            endif
          else
            if (n.eq.1) then
              rwork(1)=cabs1(b(j,1))+cabs1(d(1))*cabs1(x(j,1))
            else
              rwork(1)=cabs1(b(j,1))+cabs1(d(1))*cabs1(x(j,1))
     &          +cabs1(dl(1))*cabs1(x(j,2))
              do i=2,n - 1
                rwork(i)=cabs1(b(j,i))+cabs1(du(i-1))*cabs1(x(j,i-1))
     &            +cabs1(d(i))*cabs1(x(j,i))
     &            +cabs1(dl(i))*cabs1(x(j,i+1))
              enddo
              rwork(n)=cabs1(b(j,n))+cabs1(du(n-1))*cabs1(x(j,n-1))
     &          +cabs1(d(n))*cabs1(x(j,n))
            endif
          endif
          s=zero
          do i=1,n
            if (rwork(i).gt.safe2) then
              s=max(s,cabs1(work(i))/rwork(i))
            else
              s=max(s,(cabs1(work(i))+safe1)/(rwork(i)+safe1))
            endif
          enddo
          berr(j)=s
          if (berr(j).gt.eps .and. two*berr(j).le.lstres .and.
     &    count.le.itmax) then
            if (lsame(trans,'C')) then
              do i=1,n
                work(i)=conjg(work(i))
              enddo
            endif
            call cgttrs(ltrans,n,1,dlf,df,duf,du2,ipiv,work,n,info)
            if (lsame(trans,'C')) then
              do i=1,n
                work(i)=conjg(work(i))
              enddo
            endif
            call zaxpy(n,cmplx(one),work,1,x(j,1),ldx)
            lstres=berr(j)
            count=count+1
            goto 20
          endif
        do i=1,n
          if (rwork(i).gt.safe2) then
            rwork(i)=cabs1(work(i))+nz*eps*rwork(i)
          else
            rwork(i)=cabs1(work(i))+nz*eps*rwork(i)+safe1
          endif
        enddo
        kase=0
   70   continue
          call zlacn2(n,work(n+1),work,ferr(j),kase,isave)
          if (kase.ne.0) then
            if (kase.eq.1) then
              call cgttrs(transt,n,1,dlf,df,duf,du2,ipiv,work,n,info)
              do i=1,n
                work(i)=rwork(i)*work(i)
              enddo
            else
              do i=1,n
                work(i)=rwork(i)*work(i)
              enddo
              call cgttrs(transn,n,1,dlf,df,duf,du2,ipiv,work,n,info)
            endif
            goto 70
          endif
        lstres=zero
        do i=1,n
          lstres=max(lstres,cabs1(x(j,i)))
        enddo
        if (lstres.ne.zero) ferr(j)=ferr(j)/lstres
      enddo
      return
      end
