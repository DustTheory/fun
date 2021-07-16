c modified from zgtrfsr to avoid pivoting
      subroutine zgtrfsrnp(trans,n,nrhs,dl,d,du,dlf,df,duf,b,ldb,
     &  x,ldx,ferr,berr,work,rwork,info)
      character trans
      integer info,ldb,ldx,n,nrhs
      double precision berr(*),ferr(*),rwork(*)
      complex*16 b(ldb,*),d(*),df(*),dl(*),dlf(*),du(*),duf(*),
     &  work(*),x(ldx,*)
      integer itmax
      parameter (itmax=5)
      double precision zero,one
      parameter (zero=0.0d+0,one=1.0d+0)
      double precision two
      parameter (two=2.0d+0)
      double precision three
      parameter (three=3.0d+0)
      logical notran
      character ltrans,transn,transt
      integer count,i,j,kase,nz
      double precision eps,lstres,s,safe1,safe2,safmin
      complex*16 zdum
      integer isave(3)
      external xerbla,zaxpy,zcopy,zgttrs,zlacn2,zlagtm
      intrinsic abs,dble,dcmplx,dimag,max
      logical lsame
      double precision dlamch
      external lsame,dlamch
      double precision cabs1
      cabs1(zdum)=abs(dble(zdum))+abs(dimag(zdum))

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
        call xerbla('zgtrfsrnp',-info)
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
            enddo
          endif
c         call zlagtm(ltrans,n,1,-one,dl,d,du,x(j,1),ldx,one,work,n)
          if (notran) then
            call zgtmv(n,dcmplx(-one),du,d,dl,x(j,1),ldx,dcmplx(one),
     &        work,1)
          else
            call zgtmv(n,dcmplx(-one),dl,d,du,x(j,1),ldx,dcmplx(one),
     &        work,1)
          endif
          if (lsame(trans,'C')) then
            do i=1,n
              work(i)=conjg(work(i))
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
            call zgttrsnp(ltrans,n,1,dlf,df,duf,work,n,info)
            if (lsame(trans,'C')) then
              do i=1,n
                work(i)=conjg(work(i))
              enddo
            endif
            call zaxpy(n,dcmplx(one),work,1,x(j,1),ldx)
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
              call zgttrsnp(transt,n,1,dlf,df,duf,work,n,info)
              do i=1,n
                work(i)=rwork(i)*work(i)
              enddo
            else
              do i=1,n
                work(i)=rwork(i)*work(i)
              enddo
              call zgttrsnp(transn,n,1,dlf,df,duf,work,n,info)
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
