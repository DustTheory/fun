      function riemann(a,b,n,f)
      double precision riemann
      double precision a,b
      integer n
      double precision f
      external f
      double precision dx,x
      integer i
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      riemann=0.d0
      dx=(b-a)/dble(n)
      x=a
      do i=1,n
        x=x+dx
        riemann=riemann+f(x)
      enddo
      riemann=riemann*dx
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

c***********************************************************************

      function midpoint(a,b,n,f)
      double precision midpoint
      double precision a,b
      integer n
      double precision f
      external f
      double precision dx,x
      integer i
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      midpoint=0.d0
      dx=(b-a)/dble(n)
      x=a+0.5*dx
      do i=1,n
        midpoint=midpoint+f(x)
        x=x+dx
      enddo
      midpoint=midpoint*dx
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

c***********************************************************************

      function trapezoidal(a,b,n,f)
      double precision trapezoidal
      double precision a,b
      integer n
      double precision f
      external f
      double precision dx,x
      integer i
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      trapezoidal=0.d0
      dx=(b-a)/dble(n)
      x=a
      do i=1,n-1
        x=x+dx
        trapezoidal=trapezoidal+f(x)
      enddo
      trapezoidal=(trapezoidal+(f(a)+f(b))*0.5d0)*dx
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end
