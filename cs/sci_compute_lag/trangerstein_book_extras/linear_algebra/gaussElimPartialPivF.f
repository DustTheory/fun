      program main
      integer i,j,n,info
      parameter (n=10)
      integer ipiv(n)
      double precision A(n,n),LU(n,n),b(n),x(n)
      do j=1,n
        do i=1,n
          A(i,j)=rand()
          LU(i,j)=A(i,j)
        enddo
      enddo
      call dgetrf(n,n,LU,n,ipiv,info)
      print *, "info = ",info
      print *, "ipiv = ",(ipiv(i),i=1,n-1)

      do j=1,n
        b(j)=rand()
      enddo
      call dcopy(n,b,1,x,1)
      call dgetrs('N',n,1,LU,n,ipiv,x,n,info)
      print *, "x = ",(x(j),j=1,n)

      call dgemv('N',n,n,-1.d0,A,n,x,1,1.d0,b,1)
      print *, "error = ",(b(j),j=1,n)
      return
      end
