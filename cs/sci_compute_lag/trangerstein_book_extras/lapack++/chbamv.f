c modifed from chbmv to work with absolute values
      subroutine chbamv(uplo,n,k,alpha,a,lda,x,incx,beta,y,incy)
      real alpha,beta
      integer incx,incy,k,lda,n
      character uplo
      complex a(lda,*),x(*)
      real y(*)
      complex one
      parameter (one=(1.0,0.0))
      complex zero
      parameter (zero=(0.0,0.0))
      real temp1,temp2
      integer i,info,ix,iy,j,jx,jy,kplus1,kx,ky,l
      logical lsame
      external lsame
      external xerbla
      intrinsic max,min

      info=0
      if (.not.lsame(uplo,'U') .and. .not.lsame(uplo,'L')) then
        info=1
      else if (n.lt.0) then
        info=2
      else if (k.lt.0) then
        info=3
      else if (lda.lt. (k+1)) then
        info=6
      else if (incx.eq.0) then
        info=8
      else if (incy.eq.0) then
        info=11
      endif
      if (info.ne.0) then
        call xerbla('chbamv ',info)
        return
      endif
      if ((n.eq.0) .or. ((alpha.eq.zero).and. (beta.eq.one))) return
      if (incx.gt.0) then
        kx=1
      else
        kx=1-(n-1)*incx
      endif
      if (incy.gt.0) then
        ky=1
      else
        ky=1-(n-1)*incy
      endif
      if (beta.ne.one) then
        if (incy.eq.1) then
          if (beta.eq.zero) then
            do i=1,n
              y(i)=zero
            enddo
          else
            do i=1,n
              y(i)=beta*y(i)
            enddo
          endif
        else
          iy=ky
          if (beta.eq.zero) then
            do i=1,n
              y(iy)=zero
              iy=iy+incy
            enddo
          else
            do i=1,n
              y(iy)=beta*y(iy)
              iy=iy+incy
            enddo
          endif
        endif
      endif
      if (alpha.eq.zero) return
      if (lsame(uplo,'U')) then
        kplus1=k+1
        if ((incx.eq.1) .and. (incy.eq.1)) then
          do j=1,n
            temp1=alpha*abs(x(j))
            temp2=zero
            l=kplus1-j
            do i=max(1,j-k),j-1
              y(i)=y(i)+temp1*abs(a(l+i,j))
              temp2=temp2+abs(a(l+i,j))*abs(x(i))
            enddo
            y(j)=y(j)+temp1*abs(a(kplus1,j))+alpha*temp2
          enddo
        else
          jx=kx
          jy=ky
          do j=1,n
            temp1=alpha*abs(x(jx))
            temp2=zero
            ix=kx
            iy=ky
            l=kplus1-j
            do i=max(1,j-k),j-1
              y(iy)=y(iy)+temp1*abs(a(l+i,j))
              temp2=temp2+abs(a(l+i,j))*abs(x(ix))
              ix=ix+incx
              iy=iy+incy
            enddo
            y(jy)=y(jy)+temp1*abs(a(kplus1,j))+alpha*temp2
            jx=jx+incx
            jy=jy+incy
            if (j.gt.k) then
              kx=kx+incx
              ky=ky+incy
            endif
          enddo
        endif
      else
        if ((incx.eq.1) .and. (incy.eq.1)) then
          do j=1,n
            temp1=alpha*abs(x(j))
            temp2=zero
            y(j)=y(j)+temp1*abs(a(1,j))
            l=1-j
            do i=j+1,min(n,j+k)
              y(i)=y(i)+temp1*abs(a(l+i,j))
              temp2=temp2+abs(a(l+i,j))*abs(x(i))
            enddo
            y(j)=y(j)+alpha*temp2
          enddo
        else
          jx=kx
          jy=ky
          do j=1,n
            temp1=alpha*abs(x(jx))
            temp2=zero
            y(jy)=y(jy)+temp1*abs(a(1,j))
            l=1-j
            ix=jx
            iy=jy
            do i=j+1,min(n,j+k)
              ix=ix+incx
              iy=iy+incy
              y(iy)=y(iy)+temp1*abs(a(l+i,j))
              temp2=temp2+abs(a(l+i,j))*abs(x(ix))
            enddo
            y(jy)=y(jy)+alpha*temp2
            jx=jx+incx
            jy=jy+incy
          enddo
        endif
      endif
      return
      end
