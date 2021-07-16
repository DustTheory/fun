      subroutine gauss_elim(n,A)
      integer j,k,kp1,n
      real A(n,n)
      do k=1,n
        if (abs(A(k,k)).le.0.) return
        kp1=k+1
        call sscal_simple(n-k,1./A(k,k),A(kp1,k))
        do j=kp1,n
          call saxpy_simple(n-k,-A(k,j),A(kp1,k),A(kp1,j))
        enddo
      enddo
      return
      end
