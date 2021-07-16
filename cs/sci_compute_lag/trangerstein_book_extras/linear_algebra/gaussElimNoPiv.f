c     triangular systems overwrite A:
c       left-triangular coefficients stored in A below the diagonal
c       right-triangular coefficients stored on and above the diagonal
      subroutine gaussElimNoPiv(n,A,lda)
      integer i,j,k,kp1,lda,n
      real A(lda,*),diag
      do k=1,n
        diag=A(k,k)
        if (abs(diag).le.0.) return
        kp1=k+1
        diag=1./diag
        do i=kp1,n
          A(i,k)=A(i,k)*diag
        enddo
        do j=kp1,n
          do i=kp1,n
            A(i,j)=A(i,j)-A(i,k)*A(k,j)
          enddo
        enddo
      enddo
      return
      end

      subroutine gaussElimNoPivWithBlas(n,A,lda)
      integer i,j,k,kp1,lda,n
      real A(lda,*),diag
      do k=1,n
        diag=A(k,k)
        if (abs(diag).le.0.) return
        kp1=k+1
        call sscal(n-k,1./diag,A(kp1,k),1)
        call sger(n-k,n-k,-1.,A(kp1,k),1,A(k,kp1),lda,A(kp1,kp1),lda)
      enddo
      return
      end
