c modified from dptrfs to improve X A=B <==> A X^T=B^T
      subroutine dptrfsr(n,nrhs,d,e,df,ef,b,ldb,x,ldx,ferr,berr,
     &  work,info)
      integer info,ldb,ldx,n,nrhs
      double precision b(ldb,*),berr(*),d(*),df(*),e(*),ef(*),ferr(*),
     &  work(*),x(ldx,*)
      integer itmax
      parameter (itmax=5)
      double precision zero
      parameter (zero=0.0d+0)
      double precision one
      parameter (one=1.0d+0)
      double precision two
      parameter (two=2.0d+0)
      double precision three
      parameter (three=3.0d+0)
      integer count,i,ix,j,nz
      double precision bi,cx,dx,eps,ex,lstres,s,safe1,safe2,safmin
      external daxpy,dpttrs,xerbla
      intrinsic abs,max
      integer idamax
      double precision dlamch
      external idamax,dlamch

      info=0
      if (n.lt.0) then
        info=-1
      else if (nrhs.lt.0) then
        info=-2
      else if (ldb.lt.max(1,nrhs)) then
        info=-8
      else if (ldx.lt.max(1,nrhs)) then
        info=-10
      endif
      if (info.ne.0) then
        call xerbla('dptrfsr',-info)
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
      eps=dlamch('Epsilon')
      safmin=dlamch('Safe minimum')
      safe1=nz*safmin
      safe2=safe1/eps
      do j=1,nrhs
        count=1
        lstres=three
   20   continue
        if (n.eq.1) then
          bi=b(j,1)
          dx=d(1)*x(j,1)
          work(n+1)=bi-dx
          work(1)=abs(bi)+abs(dx)
        else
          bi=b(j,1)
          dx=d(1)*x(j,1)
          ex=e(1)*x(j,2)
          work(n+1)=bi-dx-ex
          work(1)=abs(bi)+abs(dx)+abs(ex)
          do i=2,n-1
            bi=b(j,i)
            cx=e(i-1)*x(j,i-1)
            dx=d(i)*x(j,i)
            ex=e(i)*x(j,i+1)
            work(n+i)=bi-cx-dx-ex
            work(i)=abs(bi)+abs(cx)+abs(dx)+abs(ex)
          enddo
          bi=b(j,n)
          cx=e(n-1)*x(j,n-1)
          dx=d(n)*x(j,n)
          work(n+n)=bi-cx-dx
          work(n)=abs(bi)+abs(cx)+abs(dx)
        endif
        s=zero
        do i=1,n
          if (work(i).gt.safe2) then
            s=max(s,abs(work(n+i))/work(i))
          else
            s=max(s,(abs(work(n+i))+safe1)/(work(i)+safe1))
          endif
        enddo
        berr(j)=s
        if (berr(j).gt.eps .and. two*berr(j).le.lstres .and.
     &  count.le.itmax) then
          call dpttrs(n,1,df,ef,work(n+1),n,info)
          call daxpy(n,one,work(n+1),1,x(j,1),ldx)
          lstres=berr(j)
          count=count+1
          goto 20
        endif
        do i=1,n
          if (work(i).gt.safe2) then
            work(i)=abs(work(n+i))+nz*eps*work(i)
          else
            work(i)=abs(work(n+i))+nz*eps*work(i)+safe1
          endif
        enddo
        ix=idamax(n,work,1)
        ferr(j)=work(ix)
        work(1)=one
        do i=2,n
          work(i)=one+work(i-1)*abs(ef(i-1))
        enddo
        work(n)=work(n)/df(n)
        do i=n-1,1,-1
          work(i)=work(i)/df(i)+work(i+1)*abs(ef(i))
        enddo
        ix=idamax(n,work,1)
        ferr(j)=ferr(j)*abs(work(ix))
        lstres=zero
        do i=1,n
            lstres=max(lstres,abs(x(j,i)))
        enddo
        if (lstres.ne.zero) ferr(j)=ferr(j)/lstres
      enddo
      return
      end
