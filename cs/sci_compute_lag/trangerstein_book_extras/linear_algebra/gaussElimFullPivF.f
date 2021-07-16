      program main
      integer i,j,n,info
      parameter (n=10)
      integer ipiv(n),jpiv(n)
      double precision A(n,n),LU(n,n),b(n),x(n),s
      do j=1,n
        do i=1,n
          A(i,j)=rand()
          LU(i,j)=A(i,j)
        enddo
      enddo
      call dgetc2(n,LU,n,ipiv,jpiv,info)
      print *, "info = ",info
      print *, "ipiv = ",(ipiv(i),i=1,n-1)
      print *, "jpiv = ",(jpiv(j),j=1,n-1)

      do j=1,n
        b(j)=rand()
      enddo
      call dcopy(n,b,1,x,1)
      call dgesc2(n,LU,n,x,ipiv,jpiv,s)
      print *, "scale = ",s
      print *, "x = ",(x(j),j=1,n)

      call dgemv('N',n,n,-s,A,n,x,1,1.d0,b,1)
      print *, "error = ",(b(j),j=1,n)
      return
      end
