      subroutine divdif(n,x,y, difs)
      integer n
      double precision x(0:n),y(0:n)
      double precision difs(0:n)
      integer i,j

      do j=0,n
        difs(j)=y(j)
c       print *, "difs[",j,"] = ",difs(j)
      enddo
c     print *, " "
      do j=1,n
        do i=n,j,-1
          difs(i)=(difs(i)-difs(i-1))/(x(i)-x(i-j))
c         print *, "difs[",i,"] = ",difs(i)
        enddo
c       print *, " "
      enddo

      return
      end

c***********************************************************************

      function newtonpoly(n, difs,t,x)
      double precision newtonpoly
      integer n
      double precision difs(0:n),t,x(0:n)
      integer j

      newtonpoly=difs(n)
      do j=n-1,0,-1
        newtonpoly=difs(j)+newtonpoly*(t-x(j))
      enddo

      return
      end

c***********************************************************************

      subroutine addpoint(n, x0,y0, difs,x)
      integer n
      double precision x(0:n),x0,y0
      double precision difs(0:n+1)
      integer j
      double precision temp1,temp2

      n=n+1
      do j=1,n
        x(j)=x(j-1)
      enddo
      x(0)=x0

      temp1=difs(0)
      difs(0)=y0
      do j=1,n
        temp2=temp1
        temp1=difs(j)
        difs(j)=(temp2-difs(j-1))/(x(j)-x0)
      enddo

      return
      end

c***********************************************************************

      function lagrangepoly(n, t,x,y)
      double precision lagrangepoly
      integer n
      double precision t,x(0:n),y(0:n)
      integer i,j
      double precision prod

      lagrangepoly=0.d0
      do j=0,n
        prod=y(j)
        do i=0,j-1
          prod=prod*(t-x(i))/(x(j)-x(i))
        enddo
        do i=j+1,n
          prod=prod*(t-x(i))/(x(j)-x(i))
        enddo
        lagrangepoly=lagrangepoly+prod
      enddo

      return
      end
