      program main
      integer i,j,m
      parameter (m=10)
      real L(m,m),b(m),y(m),Ly(m)
      do j=1,m
        L(j,j)=1.
        b(j)=rand()
        do i=j+1,m
          L(i,j)=rand()
        enddo
      enddo
      
c     solve L * y = b:
      call scopy(m,b,1,y,1) ! copy b to y
      call strsv('L','N','U',m,L,m,y,1) !solve L*y = b, store soln in y
      print *, "y = ",(y(j),j=1,m)

c     compute residual = b - L * y
      call scopy(m,y,1,Ly,1) ! copy y to Ly
      call strmv('L','N','U',m,L,m,Ly,1) ! multiply Ly = L * y
      call saxpy(m,-1.,Ly,1,b,1) ! resid = -Ly + b, store resid in b
      print *, "error = ",(b(i),i=1,m)

      stop
      end
