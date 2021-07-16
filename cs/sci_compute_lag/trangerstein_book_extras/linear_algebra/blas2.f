c x:= A * x or A^T * x, A unit/non-unit upper/lower triangular
c assume incx=1
      subroutine strmv_simple(uplo,trans,diag,n,A,lda,x)
      character diag,trans,uplo
      integer lda,n
      real A(lda,*),x(*)

      integer j
      real temp

      if (trans.eq.'N') then
        if (uplo.eq.'U') then
          do j=1,n
            call saxpy_simple(j-1,x(j),A(1,j),x)
            if (diag.eq.'N') x(j)=x(j)*A(j,j)
          enddo
        else
          do j=n,1,-1
            call saxpy_simple(n-j,x(j),A(j+1,j),x(j+1))
            if (diag.eq.'N') x(j)=x(j)*A(j,j)
          enddo
        endif
      else
        if (uplo.eq.'U') then
          do j=n,1,-1
            temp=x(j)
            if (diag.eq.'N') temp=temp*A(j,j)
            x(j)=temp+sdot_simple(j-1,A(1,j),x)
          enddo
        else
          do j=1,n
            temp=x(j)
            if (diag.eq.'N') temp=temp*A(j,j)
            x(j)=temp+sdot_simple(n-j,A(j+1,j),x(j+1))
          enddo
        endif
      endif
      return
      end

c x:= A^{-1} * x or A^{-T} * x, A unit/non-unit upper/lower triangular
c assume incx=1
      subroutine strsv_simple(uplo,trans,diag,n,A,lda,x)
      character diag,trans,uplo
      integer lda,n
      real A(lda,*),x(*)

      integer j
      real temp

      if (trans.eq.'N') then
        if (uplo.eq.'U') then
          do j=n,1,-1
            if (diag.eq.'N') x(j)=x(j)/A(j,j)
            call saxpy_simple(j-1,-x(j),A(1,j),x)
          enddo
        else
          do j=1,n
            if (diag.eq.'N') x(j)=x(j)/A(j,j)
            call saxpy_simple(n-j,-x(j),A(j+1,j),x(j+1))
          enddo
        endif
      else
        if (uplo.eq.'U') then
          do j=1,n
            x(j)=x(j)-sdot_simple(j-1,A(1,j),x)
            if (diag.eq.'N') x(j)=x(j)/A(j,j)
          enddo
        else
          do j=n,1,-1
            x(j)=x(j)-sdot_simple(n-j,A(j+1,j),x(j+1))
            if (diag.eq.'N') x(j)=x(j)/A(j,j)
          enddo
        endif
      endif
      return
      end

c y := A * x * alpha + y * beta or A^T * x * alpha + y * beta
c assume incx = 1 = incy
      subroutine sgemv_simple(trans,m,n,alpha,A,lda,x,beta,y)
      character trans
      integer lda,m,n
      real A(lda,*),alpha,beta,x(*),y(*)

      integer j
      real temp

      call sscal_simple(m,beta,y)
      if (trans.eq.'N') then
        do j=1,n
          call saxpy_simple(m,alpha*x(j),A(1,j),y)
        enddo
      else
        do j=1,n
          y(j)=y(j)+alpha*sdot_simple(m,A(1,j),x)
        enddo
      endif
      return
      end

c A := A + x * alpha * y^T
c assume incx = 1 = incy
      subroutine sger_simple(m,n,alpha,x,y,A,lda)
      integer lda,m,n
      real A(lda,*),x(*),y(*)

      integer j

      do j=1,n
        call saxpy_simple(m,alpha*y(j),x,A(1,j))
      enddo
      return
      end
