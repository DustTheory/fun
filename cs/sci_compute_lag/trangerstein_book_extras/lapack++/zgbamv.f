c modifed from zgbmv to use absolute values
      subroutine zgbamv(trans,m,n,kl,ku,alpha,a,lda,x,incx,beta,y,incy)
      complex*16 alpha,beta
      integer incx,incy,kl,ku,lda,m,n
      character trans
      complex*16 a(lda,*),x(*)
      double precision y(*)
      complex*16 one
      parameter (one=(1.0d+0,0.0d+0))
      complex*16 zero
      parameter (zero=(0.0d+0,0.0d+0))
      complex*16 temp
      integer i,info,ix,iy,j,jx,jy,k,kup1,kx,ky,lenx,leny
      logical noconj
      logical lsame
      external lsame
      external xerbla
      intrinsic dconjg,max,min

      info=0
      if (.not.lsame(trans,'N') .and. .not.lsame(trans,'T') .and.
     +.not.lsame(trans,'C')) then
        info=1
      else if (m.lt.0) then
        info=2
      else if (n.lt.0) then
        info=3
      else if (kl.lt.0) then
        info=4
      else if (ku.lt.0) then
        info=5
      else if (lda.lt.(kl+ku+1)) then
        info=8
      else if (incx.eq.0) then
        info=10
      else if (incy.eq.0) then
        info=13
      endif
      if (info.ne.0) then
        call xerbla('zgbamv ',info)
        return
      endif
      if ((m.eq.0) .or. (n.eq.0) .or.
     &  ((alpha.eq.zero).and. (beta.eq.one))) return
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
      if (beta.ne.one) then
        if (incy.eq.1) then
          if (beta.eq.zero) then
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
          if (beta.eq.zero) then
            do i=1,leny
              y(iy)=zero
              iy=iy+incy
            enddo
          else
            do i=1,leny
              y(iy)=beta*y(iy)
              iy=iy+incy
            enddo
          endif
        endif
      endif
      if (alpha.eq.zero) return
      kup1=ku+1
      if (lsame(trans,'n')) then
        jx=kx
        if (incy.eq.1) then
          do j=1,n
            if (x(jx).ne.zero) then
              temp=alpha*abs(x(jx))
              k=kup1-j
              do i=max(1,j-ku),min(m,j+kl)
                y(i)=y(i)+temp*abs(a(k+i,j))
              enddo
            endif
            jx=jx+incx
          enddo
        else
          do j=1,n
            if (x(jx).ne.zero) then
              temp=alpha*abs(x(jx))
              iy=ky
              k=kup1-j
              do i=max(1,j-ku),min(m,j+kl)
                y(iy)=y(iy)+temp*abs(a(k+i,j))
                iy=iy+incy
              enddo
            endif
            jx=jx+incx
            if (j.gt.ku) ky=ky+incy
          enddo
        endif
      else
        jy=ky
        if (incx.eq.1) then
          do j=1,n
            temp=zero
            k=kup1-j
            if (noconj) then
              do i=max(1,j-ku),min(m,j+kl)
                  temp=temp+abs(a(k+i,j)*x(i))
              enddo
            else
              do i=max(1,j-ku),min(m,j+kl)
                temp=temp+abs(dconjg(a(k+i,j))*x(i))
              enddo
            endif
            y(jy)=y(jy)+alpha*temp
            jy=jy+incy
          enddo
        else
          do j=1,n
            temp=zero
            ix=kx
            k=kup1-j
            if (noconj) then
              do i=max(1,j-ku),min(m,j+kl)
                temp=temp+abs(a(k+i,j)*x(ix))
                ix=ix+incx
              enddo
            else
              do i=max(1,j-ku),min(m,j+kl)
                temp=temp+abs(dconjg(a(k+i,j))*x(ix))
                ix=ix+incx
              enddo
            endif
            y(jy)=y(jy)+alpha*temp
            jy=jy+incy
            if (j.gt.ku) kx=kx+incx
          enddo
        endif
      endif
      return
      end
