      program main
      integer i,j,n,info
      character equed
      double precision rcond,ferr(1),berr(1)
      parameter (n=1024)
      integer ipiv(n),iwork(n)
      double precision A(n,n),LU(n,n),b(n),x(n)
      double precision r(n),c(n),work(4*n)
      do j=1,n
        do i=1,n
          A(i,j)=rand()
        enddo
      enddo
      do i=1,n
        b(i)=rand()
      enddo
      
      call dgesvx('E','N',n,1,A,n,LU,n,ipiv,equed,r,c,b,n,x,n,
     &  rcond,ferr,berr,work,iwork,info)
      print *, "scale and improve: info,equed,rcond,ferr,berr = ",
     &  info,equed,rcond,ferr(1),berr(1)

      return
      end
