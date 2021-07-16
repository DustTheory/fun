      subroutine cgtmv(n,alpha,L,D,U,x,incx,beta,y,incy)
c     dgtmv performs y := A * x * alpha + y * beta
c     where alpha and beta are scalars, x and y are vectors and
c     A is a tridiagonal matrix with
c     subdiagonal L, diagonal D and superdiagonal U

c     n is the order of A; n >= 0
c     incx is the stride for x
c     incy is the stride for y

      complex alpha,beta
      integer incx,incy,n
      complex L(*),D(*),U(*),x(*),y(*)

      complex temp
      integer i,info,iy,j,jx,ky
      complex one,zero
      parameter (one=(1.,0.),zero=(0.,0.))

      info=0
      if (n.lt.1) then
        info=1
      else if (incx.le.0) then
        info=7
      else if (incy.le.0) then
        info=10
      endif
      if (info.ne.0) then
        call xerbla('cgtmv',info)
        return
      endif

      if ((abs(alpha).eq.zero) .and. (beta.eq.one)) return;
      if (n.eq.1) then
        y(1)=D(1)*x(1)*alpha+y(1)*beta
        return
      endif
      if (n.eq.2) then
        y(1     )=y(1     )*beta+(D(1)*x(1)+U(1)*x(1+incx))*alpha
        y(1+incy)=y(1+incy)*beta+(L(1)*x(1)+D(2)*x(1+incx))*alpha
        return
      endif

      ky=1

      if (beta.ne.one) then
        if (incy.eq.1) then
          if (abs(beta).eq.zero) then
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
          if (abs(beta).eq.zero) then
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
      if (abs(alpha).eq.zero) return

      jx=1
      if (incy.eq.1) then
        if (abs(x(jx)).ne.zero) then
          temp=x(jx)*alpha
          y(1)=y(1)+D(1)*temp
          y(2)=y(2)+L(1)*temp
        endif
        jx=jx+incx
        do j=2,n-1
          if (abs(x(jx)).ne.zero) then
            temp=x(jx)*alpha
            y(j-1)=y(j-1)+U(j-1)*temp
            y(j)=y(j)+D(j)*temp
            y(j+1)=y(j+1)+L(j)*temp
          endif
          jx=jx+incx
        enddo
        if (abs(x(jx)).ne.zero) then
          temp=x(jx)*alpha
          y(n-1)=y(n-1)+U(n-1)*temp
          y(n)=y(n)+D(n)*temp
        endif
      else
        if (abs(x(jx)).ne.zero) then
          temp=x(jx)*alpha
          iy=ky
          y(iy)=y(iy)+D(1)*temp
          iy=iy+incy
          y(iy)=y(iy)+L(1)*temp
        endif
        jx=jx+incx
        do j=2,n-1
          if (abs(x(jx)).ne.zero) then
            iy=ky
            temp=x(jx)*alpha
            y(iy)=y(iy)+U(j-1)*temp
            iy=iy+incy
            y(iy)=y(iy)+D(j)*temp
            iy=iy+incy
            y(iy)=y(iy)+L(j)*temp
          endif
          jx=jx+incx
          ky=ky+incy
        enddo
        if (abs(x(jx)).ne.zero) then
          temp=x(jx)*alpha
          iy=ky
          y(iy)=y(iy)+U(n-1)*temp
          iy=iy+incy
          y(iy)=y(iy)+D(n)*temp
        endif
      endif

      return
      end
