      program main
      integer :: i,j,n
      real, allocatable :: R(:,:),y(:),x(:),Rx(:)

      n=10
      allocate(R(n,n))
      do j=1,n
        do i=1,j
          R(i,j)=rand() 
        enddo
      enddo

      allocate(y(n))
      call random_number(y)

      allocate(x(n));
      x=y
      call strsv('U','N','N',n,R,n,x,1)
      print *, "x = ",(x(j),j=1,n)

      allocate(Rx(n))
      Rx=x
      call strmv('U','N','N',n,R,n,Rx,1)
      y=y-Rx
      print *, "error = ",(y(i),i=1,n)

      deallocate(Rx)
      deallocate(x)
      deallocate(y)
      deallocate(R)

      stop
      end
