c modified from dgtrfs to multiply by matrix on the right and avoid pivoting
      subroutine dgtrfsrnp(trans,n,nrhs,dl,d,du,dlf,df,duf,b,ldb,x,ldx,
     &  ferr,berr,work,iwork,info)
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
      integer count, i,j,kase,nz
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
     &.not.lsame(trans,'C')) then
        info=-1
      else if (n.lt.0) then
        info=-2
      else if (nrhs.lt.0) then
        info=-3
      else if (ldb.lt.max(1,nrhs)) then
        info=-11
      else if (ldx.lt.max(1,nrhs)) then
        info=-13
      endif
      if (info.ne.0) then
        call xerbla('dgtrfsrnp',-info)
        return
      endif
      if (n.eq.0 .or. nrhs.eq.0) then
        do j = 1, nrhs
          ferr(j)=zero
          berr(j)=zero
        enddo
        return
      end if
      if (.not.notran) then
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
          if (.not.notran) then ! X A^T = B <==> A X^T = B^T
            if (n.eq.1) then
              work(1)=abs(b(j,1))+abs(d(1)*x(j,1))
              work(n+1)=b(j,1)-d(1)*x(j,1)
            else
              work(1)=abs(b(j,1))+abs(d(1)*x(j,1))+abs(du(1)*x(j,2))
              work(n+1)=b(j,1)+d(1)*x(j,1)+du(1)*x(j,2)
              do i=2,n-1
                work(i)=abs(b(j,i))+abs(dl(i-1)*x(j,i-1))
     &                 +abs(d(i)*x(j,i))+abs(du(i)*x(j,i+1))
                work(n+i)=b(j,i)+dl(i-1)*x(j,i-1)
     &                   +d(i)*x(j,i)+du(i)*x(j,i+1)
              enddo
              work(n)=abs(b(j,n))+abs(dl(n-1)*x(j,n-1))+abs(d(n)*x(j,n))
              work(n+n)=b(j,n)+dl(n-1)*x(j,n-1)+(n)*x(j,n)
           endif
         else ! X A = B <==> A^T X^T = B^T
           if (n.eq.1) then
             work(1)=abs(b(j,1))+abs(d(1)*x(j,1))
             work(n+1)=b(j,1)+d(1)*x(j,1)
           else
             work(1)=abs(b(j,1))+abs(d(1)*x(j,1))+abs(dl(1)*x(j,2))
             work(n+1)=b(j,1)+d(1)*x(j,1)+dl(1)*x(j,2)
             do i=2,n-1
               work(i)=abs(b(j,i))+abs(du(i-1)*x(j,i-1))
     &                +abs(d(i)*x(j,i))+abs(dl(i)*x(j,i+1))
               work(n+i)=b(j,i)+du(i-1)*x(j,i-1)
     &                  +d(i)*x(j,i)+dl(i)*x(j,i+1)
             enddo
             work(n)=abs(b(j,n))+abs(du(n-1)*x(j,n-1))+abs(d(n)*x(j,n))
             work(n+n)=b(j,n)+du(n-1)*x(j,n-1)+d(n)*x(j,n)
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
             endif
             goto 70
           endif
         lstres=zero
         do i= 1,n
            lstres= max(lstres,abs(x(j,i)))
         enddo
         if (lstres.ne.zero) ferr(j)=ferr(j)/lstres
      enddo
      return
      end
