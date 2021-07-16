c modified from cptrfs for X A = B <==> A X^H = B^H
      subroutine cptrfsr(uplo,n,nrhs,d,e,df,ef,b,ldb,x,ldx,ferr,berr,
     &  work,rwork,info)
      character uplo
      integer info,ldb,ldx,n,nrhs
      real berr(*),d(*),df(*),ferr(*),rwork(*)
      complex b(ldb,*),e(*),ef(*),work(*),x(ldx,*)
      integer itmax
      parameter (itmax=5)
      real zero
      parameter (zero=0.0)
      real one
      parameter (one=1.0)
      real two
      parameter (two=2.0)
      real three
      parameter (three=3.0)
      logical upper
      integer count,i,ix,j,nz
      real eps,lstres,s,safe1,safe2,safmin
      complex bi,cx,dx,ex,zdum
      logical lsame
      integer isamax
      real slamch
      external lsame,isamax,slamch
      external xerbla,caxpy,cpttrs
      real cabs1
      cabs1(zdum)=abs(real(zdum))+abs(aimag(zdum))

      info=0
      upper=lsame(uplo,'U')
      if (.not.upper .and. .not.lsame(uplo,'L')) then
        info=-1
      else if (n.lt.0) then
        info=-2
      else if (nrhs.lt.0) then
        info=-3
      else if (ldb.lt.max(1,nrhs)) then
        info=-9
      else if (ldx.lt.max(1,nrhs)) then
        info=-11
      endif
      if (info.ne.0) then
        call xerbla('cptrfsr',-info)
        return
      endif
      if (n.eq.0 .or. nrhs.eq.0) then
        do j=1,nrhs
          ferr(j)=zero
          berr(j)=zero
        enddo
        return
      endif
      nz=4
      eps=slamch('Epsilon')
      safmin=slamch('Safe minimum')
      safe1=nz*safmin
      safe2=safe1 / eps
      do j=1,nrhs
        count=1
        lstres=three
   20   continue
        if (upper) then
          if (n.eq.1) then
            bi=conjg(b(j,1))
            dx=d(1)*conjg(x(j,1))
            work(1)=bi - dx
            rwork(1)=cabs1(bi)+cabs1(dx)
          else
            bi=conjg(b(j,1))
            dx=d(1)*conjg(x(j,1))
            ex=e(1)*conjg(x(j,2))
            work(1)=bi - dx - ex
            rwork(1)=cabs1(bi)+cabs1(dx)+cabs1(e(1))*cabs1(x(j,2))
            do i=2,n - 1
              bi=conjg(b(j,i))
              cx=conjg(e(i-1))*conjg(x(j,i-1))
              dx=d(i)*conjg(x(j,i))
              ex=e(i)*conjg(x(j,i+1))
              work(i)=bi - cx - dx - ex
              rwork(i)=cabs1(bi)+cabs1(e(i-1))*cabs1(x(j,i-1))
     &          +cabs1(dx)+cabs1(e(i))*cabs1(x(i+1,j))
            enddo
            bi=conjg(b(j,n))
            cx=conjg(e(n-1))*conjg(x(j,n-1))
            dx=d(n)*conjg(x(j,n))
            work(n)=bi - cx - dx
            rwork(n)=cabs1(bi)+cabs1(e(n-1))*cabs1(x(j,n-1))+cabs1(dx)
          endif
        else
          if (n.eq.1) then
            bi=conjg(b(j,1))
            dx=d(1)*conjg(x(j,1))
            work(1)=bi - dx
            rwork(1)=cabs1(bi)+cabs1(dx)
          else
            bi=conjg(b(j,1))
            dx=d(1)*conjg(x(j,1))
            ex=conjg(e(1))*conjg(x(j,2))
            work(1)=bi - dx - ex
            rwork(1)=cabs1(bi)+cabs1(dx)+cabs1(e(1))*cabs1(x(j,2))
            do i=2,n - 1
              bi=conjg(b(j,i))
              cx=e(i-1)*conjg(x(j,i-1))
              dx=d(i)*conjg(x(j,i))
              ex=conjg(e(i))*conjg(x(j,i+1))
              work(i)=bi - cx - dx - ex
              rwork(i)=cabs1(bi)+cabs1(e(i-1))*cabs1(x(j,i-1))
     &          +cabs1(dx)+cabs1(e(i))*cabs1(x(j,i+1))
            enddo
            bi=conjg(b(j,n))
            cx=e(n-1)*conjg(x(j,n-1))
            dx=d(n)*conjg(x(j,n))
            work(n)=bi - cx - dx
            rwork(n)=cabs1(bi)+cabs1(e(n-1))*cabs1(x(j,n-1))+cabs1(dx)
          endif
        endif
        s=zero
        do i=1,n
          if (rwork(i).gt.safe2) then
            s=max(s,cabs1(work(i)) / rwork(i))
          else
            s=max(s,(cabs1(work(i))+safe1)/(rwork(i)+safe1))
          endif
        enddo
        berr(j)=s
        if (berr(j).gt.eps .and. two*berr(j).le.lstres .and.
     &  count.le.itmax) then
          call cpttrs(uplo,n,1,df,ef,work,n,info)
c         call caxpy(n,cmplx(one),work,1,x(1,j),1)
          do i=1,n
            x(j,i)=x(j,i)+conjg(work(i))
          enddo
          lstres=berr(j)
          count=count+1
          go to 20
        endif
        do i=1,n
          if (rwork(i).gt.safe2) then
            rwork(i)=cabs1(work(i))+nz*eps*rwork(i)
          else
            rwork(i)=cabs1(work(i))+nz*eps*rwork(i)+safe1
          endif
        enddo
        ix=isamax(n,rwork,1)
        ferr(j)=rwork(ix)
        rwork(1)=one
        do i=2,n
            rwork(i)=one+rwork(i-1)*abs(ef(i-1))
        enddo
        rwork(n)=rwork(n)/df(n)
        do i=n-1,1,-1
          rwork(i)=rwork(i)/df(i)+rwork(i+1)*abs(ef(i))
        enddo
        ix=isamax(n,rwork,1)
        ferr(j)=ferr(j)*abs(rwork(ix))
        lstres=zero
        do i=1,n
          lstres=max(lstres,abs(x(j,i)))
        enddo
        if (lstres.ne.zero) ferr(j)=ferr(j)/lstres
      enddo
      return
      end
