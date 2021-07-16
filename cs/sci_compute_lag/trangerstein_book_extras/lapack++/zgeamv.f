c modified from zgemv to use absolute values
      subroutine zgeamv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
      double precision alpha,beta
      integer incx,incy,lda,m,n
      character trans
      complex*16 a(lda,*),x(*)
      double precision y(*)
      complex*16 zero
      parameter (zero=(0.0d+0,0.0d+0))
      double precision temp
      integer i,info,ix,iy,j,jx,jy,kx,ky,lenx,leny
      logical noconj
      logical lsame
      external lsame
      external xerbla
      intrinsic dconjg,max

      info=0
      if (.not.lsame(trans,'N') .and. .not.lsame(trans,'T') .and.
     &.not.lsame(trans,'C')) then
        info=1
      else if (m.lt.0) then
        info=2
      else if (n.lt.0) then
        info=3
      else if (lda.lt.max(1,m)) then
        info=6
      else if (incx.eq.0) then
        info=8
      else if (incy.eq.0) then
        info=11
      endif
      if (info.ne.0) then
        call xerbla('zgeamv ',info)
        return
      endif
      if ((m.eq.0) .or. (n.eq.0) .or.
     &((alpha.eq.0.d0).and. (beta.eq.1.d0))) return
      noconj=lsame(trans,'T')
      if (lsame(trans,'N')) then
        lenx=n
        leny=m
      else
        lenx=m
        leny=n
      endif
      if (incx.gt.0) then
        kx=1
      else
        kx=1-(lenx-1)*incx
      endif
      if (incy.gt.0) then
        ky=1
      else
        ky=1-(leny-1)*incy
      endif
      if (beta.ne.1.d0) then
        if (incy.eq.1) then
          if (beta.eq.0.d0) then
            do i=1,leny
              y(i)=zero
            enddo
          else
            do i=1,leny
              y(i)=beta*y(i)
            enddo
          endif
        else
          iy=ky
          if (beta.eq.0.d0) then
            do i=1,leny
              y(iy)=zero
              iy=iy + incy
            enddo
          else
            do i=1,leny
              y(iy)=beta*y(iy)
              iy=iy + incy
            enddo
          endif
        endif
      endif
      if (alpha.eq.0.d0) return
      if (lsame(trans,'N')) then
        jx=kx
        if (incy.eq.1) then
          do j=1,n
            if (x(jx).ne.zero) then
              temp=alpha*abs(x(jx))
              do i=1,m
                y(i)=y(i)+temp*abs(a(i,j))
              enddo
            endif
            jx=jx+incx
          enddo
        else
          do j=1,n
            if (x(jx).ne.zero) then
              temp=alpha*abs(x(jx))
              iy=ky
              do i=1,m
                y(iy)=y(iy)+temp*abs(a(i,j))
                iy=iy+incy
              enddo
            endif
            jx=jx+incx
          enddo
        endif
      else
        jy=ky
        if (incx.eq.1) then
          do j=1,n
            temp=0.d0
            if (noconj) then
              do i=1,m
                  temp=temp+abs(a(i,j)*x(i))
              enddo
            else
              do i=1,m
                temp=temp+abs(dconjg(a(i,j))*x(i))
              enddo
            endif
            y(jy)=y(jy)+alpha*temp
            jy=jy+incy
          enddo
        else
          do j=1,n
            temp=0.d0
            ix=kx
            if (noconj) then
              do i=1,m
                temp=temp+abs(a(i,j)*x(ix))
                ix=ix+incx
              enddo
            else
              do i=1,m
                temp=temp+abs(dconjg(a(i,j))*x(ix))
                ix=ix+incx
              enddo
            endif
            y(jy)=y(jy)+alpha*temp
            jy=jy+incy
          enddo
        endif
      endif
      return
      end
