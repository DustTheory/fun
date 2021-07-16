c modifed from cgtts2 to avoid pivoting
      subroutine cgtts2np(itrans,n,nrhs,dl,d,du,b,ldb)
      integer itrans,ldb,n,nrhs
      complex b(ldb,*),d(*),dl(*),du(*)
      integer i, j
      intrinsic conjg
      if (n.eq.0 .or. nrhs.eq.0) return
      if (itrans.eq.0) then
        if (nrhs.le.1) then
          j=1
   10     continue
            do i=1,n-1
              b(i+1,j)=b(i+1,j)-dl(i)*b(i,j)
            enddo
            b(n,j)=b(n,j)/d(n)
            if (n.gt.1) b(n-1,j)=(b(n-1,j)-du(n-1)*b(n,j))/d(n-1)
            do i=n-2,1,-1
              b(i,j)=(b(i,j)-du(i)*b(i+1,j))/d(i)
            enddo
            if (j.lt.nrhs) then
              j=j+1
              goto 10
            endif
        else
          do j=1,nrhs
            do i=1,n-1
              b(i+1,j)=b(i+1,j)-dl(i)*b(i,j)
            enddo
            b(n,j)=b(n,j)/d(n)
            if (n.gt.1) b(n-1,j)=(b(n-1,j)-du(n-1)*b(n,j))/d(n-1)
            do i=n-2, 1, -1
              b(i,j)=(b(i,j)-du(i)*b(i+1,j))/d(i)
            enddo
          enddo
        endif
      else if (itrans.eq.1) then
        if (nrhs.le.1) then
          j=1
   70     continue
            b(1,j)=b(1,j)/d(1)
            if (n.gt.1) b(2,j)=(b(2,j)-du(1)*b(1,j))/d(2)
            do i=3,n
              b(i,j)=(b(i,j)-du(i-1)*b(i-1,j))/d(i)
            enddo
            do i=n-1,1,-1
              b(i,j)=b(i,j)-dl(i)*b(i+1,j)
            enddo
            if ( j.lt.nrhs ) then
              j=j+1
              goto 70
            endif
        else
          do j=1,nrhs
            b(1,j)=b(1,j)/d(1)
            if (n.gt.1) b(2,j)=(b(2,j)-du(1)*b(1,j))/d(2)
            do i=3,n
              b(i,j)=(b(i,j)-du(i-1)*b(i-1,j))/d(i)
            enddo
            do i=n-1,1,-1
              b(i,j)=b(i,j)-dl(i)*b(i+1,j)
            enddo
          enddo
        endif
      else
        if (nrhs.le.1) then
          j=1
  130     continue
            b(1,j)=b(1,j)/conjg(d(1))
            if (n.gt.1) 
     &        b(2,j)=(b(2,j)-conjg(du(1))*b(1,j))/conjg(d(2))
            do i=3,n
              b(i,j)=(b(i,j)-conjg(du(i-1))*b(i-1,j))/conjg(d(i))
            enddo
            do i=n-1,1,-1
              b(i,j)=b(i,j)-conjg(dl(i))*b(i+1,j)
            enddo
            if (j.lt.nrhs) then
              j=j+1
              goto 130
            endif
        else
          do j=1,nrhs
            b(1,j)=b(1,j)/conjg(d(1))
            if (n.gt.1)
     &        b(2,j)=(b(2,j)-conjg(du(1))*b(1,j))/conjg(d(2))
            do i=3,n
              b(i,j)=(b(i,j)-conjg(du(i-1))*b(i-1,j))/conjg(d(i))
            enddo
            do i=n-1,1,-1
              b(i,j)=b(i,j)-conjg(dl(i))*b(i+1,j)
            enddo
          enddo
        endif
      endif
      end
