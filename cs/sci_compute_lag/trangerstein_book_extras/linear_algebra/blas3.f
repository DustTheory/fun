c C := op(A) * alpha * op(B) + C * beta
      subroutine sgemm_simple(transa,transb,m,n,k,alpha,A,lda,
     &  B,ldb,C,ldc)
      character trans
      integer lda,ldb,ldc,k,m,n
      real A(lda,*),alpha,beta,B(ldb,*),C(ldc,*)

      integer j,l
      real temp

      if (transb.eq.'N') then
        do j=1,n
          call sgemv_simple(transa,m,k,alpha,A,lda,B(1,j),beta,C(1,j))
        enddo
      else
        do j=1,n
          call sscal_simple(m,beta,C(1,j))
          do l=1,k
            call saxpy_simple(m,alpha*B(j,l),A(1,l),C(1,j))
          enddo
        enddo
      endif
      return
      end

c B:= op(A) * alpha * B or B * alpha * op(A)
      subroutine strmm_simple(side,uplo,transa,diag,m,n,A,lda,B,ldb)
      character diag,side,transa,uplo
      integer lda,ldb,m,n
      real A(lda,*),B(ldb,)

      integer j,k
      real temp

      if (side.eq.'L') then
        do j=1,n
          call strmv_simple(uplo,transa,diag,m,A,lda,B(1,j))
          call sscal_simple(m,alpha,B(1,j))
        enddo
      else
        if (transa.eq.'N') then
          if (uplo.eq.'U') then
            do j=n,1,-1
              temp=alpha
              if (diag.eq.'N') temp=temp*A(j,j)
              call sscal(m,temp,B(1,j))
              do k=1,j-1
                temp=alpha*A(k,j)
                call saxpy(m,temp,B(1,k),B(1,j))
              enddo
            enddo
          else
            do j=1,n
              temp=alpha
              if (diag.eq.'N') temp=temp*A(j,j)
              call sscal(m,temp,B(1,j))
              do k=j+1,n
                temp=alpha*A(k,j)
                call saxpy(m,temp,B(1,k),B(1,j))
              enddo
            enddo
          endif
        else
          if (uplo.eq.'U') then
            do k=1,n
              do j=1,k-1
                call saxpy(m,alpha*A(j,k),B(1,k),B(1,j))
              enddo
              temp=alpha
              if ((diag.eq.'N') temp=temp*A(k,k)
              call sscal(m,temp,B(1,k))
            enddo
          else
            do k=n,1,-1
              do j=k+1,n
                call saxpy(m,alpha*A(j,k),B(1,k),B(1,j))
              enddo
              temp=alpha
              if ((diag.eq.'N') temp=temp*A(k,k)
              call sscal(m,temp,B(1,k))
            enddo
          endif
        endif
      endif
      return
      end

c X:= op(A) * B *alpha or B * op(A) * alpha,
c   A unit/non-unit upper/lower triangular
      subroutine strsm_simple(side,uplo,transa,diag,m,n,alpha,A,lda,
     &  B,ldb)
      character diag,side,transa,uplo
      integer lda,ldb,m,n
      real A(lda,*),B(ldb,*)

      integer j
      real temp

      if (side.eq.'L') then
        do ( j=1,n
          call strsv_simple(uplo,transa,diag,m,A,lda,B(1,j))
        enddo
      else
        if (transa.eq.'N') then
          if (uplo.eq.'U') then
            do j=1,n
              call sscal_simple(m,alpha,B(1,j))
              do k=1,j-1
                call saxpy_simple(m,-A(k,j),B(1,k),B(1,j))
              enddo
              if (diag.eq.'N') then
                call sscal_simple(m,1./A(j,j),B(1,j))
              endif
            enddo
          else
            do j=n,1,-1
              call sscal_simple(m,alpha,B(1,j))
              do k=j+1,n
                call saxpy_simple(m,-A(k,j),B(1,k),B(1,j))
              enddo
              if (diag.eq.'N') then
                call sscal_simple(m,1./A(j,j),B(1,j))
              endif
            enddo
          endif
        else
          if (uplo.eq.'U') then
            do k=n,1,-1
              if (diag.eq.'N') then
                call sscal_simple(m,1./A(k,k),B(1,k))
              endif
              do j=1,k-1
                call saxpy_simple(m,-A(j,k),B(1,k),B(1,j))
              enddo
              call sscal_simple(m,alpha,B(1,k))
            enddo
          else
            do k=1,n
              if (diag.eq.'N') then
                call sscal_simple(m,1./A(k,k),B(1,k))
              endif
              do j=k+1,n
                call saxpy_simple(m,-A(j,k),B(1,k),B(1,j))
              enddo
              call sscal_simple(m,alpha,B(1,k))
            enddo
          endif
        endif
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
