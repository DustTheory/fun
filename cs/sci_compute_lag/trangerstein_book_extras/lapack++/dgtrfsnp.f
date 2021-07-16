c     modified from dgtrfs to avoid pivoting
      subroutine dgtrfsnp(trans,n,nrhs,dl,d,du,dlf,df,duf,b,ldb,x,ldx,
     &  ferr,berr,work,iwork,info )
      character trans
      integer info,ldb,ldx,n,nrhs
      integer iwork(*)
      double precision b(ldb,*),berr(*),d(*),df(*),dl(*),dlf(*),du(*),
     &  duf(*),ferr(*),work(*),x(ldx,*)
      integer itmax
      parameter (itmax=5)
      double precision zero,one
      parameter (zero=0.0d+0,one=1.0d+0)
      double precision two
      parameter (two=2.0d+0)
      double precision three
      parameter (three=3.0d+0)
      logical notran
      character transn,transt
      integer count,i,j,kase,nz
      double precision eps,lstres,s,safe1,safe2,safmin
      integer isave(3)
      external daxpy,dcopy,dgttrs,dlacn2,dlagtm,xerbla
      intrinsic abs,max
      logical lsame
      double precision dlamch
      external lsame,dlamch

      info=0
      notran=lsame(trans,'N')
      if (.not.notran .and. .not.lsame(trans,'T') .and.
     &.not.  lsame(trans,'C')) then
        info=-1
      else if (n.lt.0) then
        info=-2
      else if (nrhs.lt.0) then
        info=-3
      else if (ldb.lt.max(1,n)) then
        info=-11
      else if (ldx.lt.max(1,n)) then
        info=-13
      end if
      if (info.ne.0) then
        call xerbla('dgtrfsnp',-info)
        return
      end if
      if (n.eq.0 .or. nrhs.eq.0) then
        do j = 1, nrhs
          ferr(j) = zero
          berr(j) = zero
        enddo
        return
      endif
      if (notran) then
        transn='N'
        transt='T'
      else
        transn='T'
        transt='N'
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
          call dcopy(n,b(1,j),1,work(n+1),1)
          call dlagtm(trans,n,1,-one,dl,d,du,x(1,j),ldx,one,work(n+1),n)
          if (notran) then
            if (n.eq.1) then
              work(1)=abs(b(1,j))+abs(d(1)*x(1,j))
            else
              work(1)=abs(b(1,j))+abs(d(1)*x(1,j))+abs(du(1)*x(2,j))
              do i=2,n-1
                  work(i)=abs(b(i,j))+abs(dl(i-1)*x(i-1,j))
     &                   +abs(d(i)*x(i,j))+abs(du(i)*x(i+1,j))
              enddo
              work(n)=abs(b(n,j))+abs(dl(n-1)*x(n-1,j))+abs(d(n)*x(n,j))
            endif
          else
            if (n.eq.1) then
              work(1)=abs(b(1,j))+abs(d(1)*x(1,j))
            else
              work(1)=abs(b(1,j))+abs(d(1)*x(1,j))+abs(dl(1)*x(2,j))
              do i = 2, n - 1
                work(i)=abs(b(i,j))+abs(du(i-1)*x(i-1,j))
     &                 +abs(d(i)*x(i,j))+abs(dl(i)*x(i+1,j))
              enddo
              work(n)=abs(b(n,j))+abs(du(n-1)*x(n-1,j))+abs(d(n)*x(n,j))
            endif
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
     &   count.le.itmax) then
           call dgttrsnp(trans,n,1,dlf,df,duf,work(n+1),n,info)
           call daxpy(n,one,work(n+1),1,x(1,j),1)
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
         kase=0
   70    continue
           call dlacn2(n,work(2*n+1),work(n+1),iwork,ferr(j),kase,isave)
           if (kase.ne.0) then
            if (kase.eq.1) then
              call dgttrsnp(transt,n,1,dlf,df,duf,work(n+1),n,info)
              do i=1,n
                work(n+i)=work(i)*work(n+i)
              enddo
            else
              do i=1,n
                work(n+i)=work(i)*work(n+i)
              enddo
              call dgttrsnp(transn,n,1,dlf,df,duf,work(n+1),n,info)
            end if
            goto 70
         endif
         lstres=zero
         do i=1,n
           lstres=max(lstres,abs(x(i,j)))
         enddo
         if (lstres.ne.zero) ferr(j)=ferr(j)/lstres
      enddo
      return
      end
