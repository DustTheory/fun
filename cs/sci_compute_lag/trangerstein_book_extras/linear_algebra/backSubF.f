      program main
      integer i,j,n
      parameter (n=10)
      real R(n,n),y(n),x(n),Rx(n)
      do j=1,n
        y(j)=rand()
        do i=1,j
          R(i,j)=rand() 
        enddo
      enddo

c     solve R * x = y:
      call scopy(n,y,1,x,1)
      call strsv('U','N','N',n,R,n,x,1)
      print *, "x = ",(x(j),j=1,n)

c     compute residual = y - R * x
      call scopy(n,x,1,Rx,1)
      call strmv('U','N','N',n,R,n,Rx,1)
      call saxpy(n,-1.,Rx,1,y,1)
      print *, "error = ",(y(i),i=1,n)

      stop
      end
