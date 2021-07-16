      function c0_spline_eval(element,nelements,order,pt,x,xi,y,p,q)
      double precision c0_spline_eval
      integer element,nelements,order
      double precision pt
      double precision x(0:nelements),xi(0:order),y(0:order*nelements)
      double precision p(0:order),q(order)

      double precision dx,t,xii
      integer i,j,k
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      j=element*order
c     print *, "in c0_spline_eval"
c     print *, "order,element,pt",order,element,pt
c     do i=element,element+1
c       print *, "x[",i,"] = ",x(i)
c     enddo
c     do i=j,j+order
c       print *, "y[",i,"] = ",y(i)
c     enddo
c     print *, "p = "
c     call print_loc(p)
c     print *, "q = "
c     call print_loc(q)
c     call mem_check()

      dx=x(element+1)-x(element)
      xii=(pt-x(element))/dx
c     print *, "xii = ",xii
      if (xii.lt.0.d0 .or. xii.gt.1.d0) then
        print *, "in c0_spline_eval"
        print *, "pt = ",pt
        print *, "element = ",element
        print *, "x = ",x(element),x(element+1)
        call abort()
      endif
c     Neville-Aitken algorithm
      do i=0,order
        p(i)=y(j+i)
      enddo
      do j=0,order-1
        t=xii-xi(j)
        do k=j+1,order
          q(k)=((xii-xi(k))*p(j)-t*p(k))/(xi(j)-xi(k))
        enddo
        do k=j+1,order
          p(k)=q(k)
        enddo
      enddo
      c0_spline_eval=p(order)

c     call mem_check()
c     print *, "leaving c0_spline_eval"
c     call flush(6)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

      subroutine c1_quadratic_coefs(nelements, x,y,yp, z)
      integer nelements
      double precision x(0:nelements),y(0:nelements),yp(0:0)
      double precision z(0:nelements)
      integer i
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      z(0)=yp(0)
      do i=1,nelements
        z(i)=2*(y(i)-y(i-1))/(x(i)-x(i-1))-z(i-1)
      enddo
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

      function c1_quadratic_eval(element,nelements,pt,x,y,z)
      double precision c1_quadratic_eval
      integer element,nelements
      double precision pt
      double precision x(0:nelements),y(0:nelements),z(0:nelements)

      double precision dx,xi
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "in c1_quadratic_eval"
c     print *, "element,pt",element,pt
c     do i=element,element+1
c       print *, "x[",i,"] = ",x(i)
c     enddo
c     do i=element,element+1
c       print *, "y[",i,"] = ",y(i)
c     enddo
      dx=x(element+1)-x(element)
      xi=(pt-x(element))/dx
c     print *, "xi = ",xi
      if (xi.lt.0.d0 .or. xi.gt.1.d0) then
        print *, "in c1_quadratic_eval"
        print *, "pt = ",pt
        print *, "element = ",element
        print *, "x = ",x(element),x(element+1)
        call abort()
      endif
      c1_quadratic_eval=y(element)
     &   +0.5*dx*xi*(z(element)*(2.d0-xi)+z(element+1)*xi)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

      function c1_spline_eval(element,nelements,order,pt,sigma,
     &  x,xi,y,yp)
      double precision c1_spline_eval
      integer element,nelements,order
      double precision pt,sigma
      double precision x(0:nelements),xi(0:order-2),
     &  y(0:nelements*(order-2)),yp(0:nelements)

      double precision dx,onemxii,prod,prod0,prod1,xii
      integer i,j,k
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      j=element*(order-2)
c     print *, "in c1_spline_eval"
c     print *, "element,pt,sigma",element,pt,sigma
c     do i=element,element+1
c       print *, "x[",i,"] = ",x(i)
c     enddo
c     do i=j,j+order-2
c       print *, "y[",i,"] = ",y(i)
c     enddo
c     do i=element,element+1
c       print *, "yp[",i,"] = ",yp(i)
c     enddo
      dx=x(element+1)-x(element)
      xii=(pt-x(element))/dx
c     print *, "xii = ",xii

      if (xii.lt.0.d0 .or. xii.gt.1.d0) then
        print *, "in c1_spline_eval"
        print *, "pt = ",pt
        print *, "element = ",element
        print *, "x = ",x(element),x(element+1)
        call abort()
      endif
      onemxii=1.d0-xii
      prod0=1.d0
      prod1=1.d0
      do i=1,order-3
        prod0=prod0*(1.d0-xii/xi(i))
        prod1=prod1*(1.d0-onemxii/xi(i))
      enddo
c     print *, "prod0,prod1 = ",prod0,prod1
      c1_spline_eval=y(j)*(1.d0+sigma*xii)*prod0*onemxii**2
     &   +y(j+order-2)*(1.d0+sigma*onemxii)*prod1*xii**2
     &   +(yp(element)*xii*prod0*onemxii**2
     &    -yp(element+1)*onemxii*prod1*xii**2)*dx
      do i=1,order-3
        prod=((xii/xi(i))*(onemxii/(1.-xi(i))))**2
        do k=1,order-3
          if (k.ne.i) then
            prod=prod*(xii-xi(k))/(xi(i)-xi(k))
          endif
        enddo
        c1_spline_eval=c1_spline_eval+y(j+i)*prod
      enddo
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

      subroutine c2_cubic_coefs(nelements, d,e,x,y,yp, z)
      integer nelements
      double precision d(0:nelements),e(nelements),
     &  x(0:nelements),y(0:nelements),yp(0:1)
      double precision z(0:nelements)
      integer i,info
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "in c2_cubic_coefs"
c     print *, "nelements = ",nelements
c     do i=0,nelements
c       print *, "x[",i,"] = ",x(i)
c     enddo
c     do i=0,nelements
c       print *, "y[",i,"] = ",y(i)
c     enddo
c     do i=0,1
c       print *, "yp[",i,"] = ",yp(i)
c     enddo
      do i=1,nelements
        e(i)=x(i)-x(i-1)
        z(i)=(y(i)-y(i-1))/(x(i)-x(i-1))
      enddo
c     do i=1,nelements
c       print *, "dx[",i,"] = ",e(i)
c     enddo
c     do i=1,nelements
c       print *, "divdif[",i,"] = ",z(i)
c     enddo
      d(0)=2*e(1)
      z(0)=(z(1)-yp(0))*6.d0
      do i=1,nelements-1
        d(i)=2*(e(i+1)+e(i))
        z(i)=(z(i+1)-z(i))*6.d0
      enddo
      d(nelements)=2*e(nelements)
      z(nelements)=(yp(1)-z(nelements))*6.d0
c     do i=0,nelements
c       print *, "d[",i,"] = ",d(i)
c     enddo
c     do i=1,nelements
c       print *, "e[",i,"] = ",e(i)
c     enddo
c     do i=0,nelements
c       print *, "z[",i,"] = ",z(i)
c     enddo
      call dptsv(nelements+1,1,d,e,z,nelements+1,info)
c     do i=0,nelements
c       print *, "z[",i,"] = ",z(i)
c     enddo
      if (info.ne.0) call abort()
c     print *, "leaving c2_cubic_coefs"
c     call flush(6)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

      function c2_cubic_eval(element,nelements,pt, x,y,z)
      double precision c2_cubic_eval
      integer element,nelements
      double precision pt
      double precision x(0:nelements),y(0:nelements),
     &  z(0:nelements)

      double precision dx,xi
c     integer i
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "in c2_cubic_eval"
c     print *, "element,pt",element,pt
c     do i=element,element+1
c       print *, "x[",i,"] = ",x(i)
c     enddo
c     do i=element,element+1
c       print *, "y[",i,"] = ",y(i)
c     enddo
      dx=x(element+1)-x(element)
      xi=(pt-x(element))/dx
c     print *, "xi = ",xi
      if (xi.lt.0.d0 .or. xi.gt.1.d0) then
        print *, "in c2_cubic_eval"
        print *, "pt = ",pt
        print *, "element = ",element
        print *, "x = ",x(element),x(element+1)
        call abort()
      endif
      c2_cubic_eval=y(element)+xi*(y(element+1)-y(element)
     & -(1.d0-xi)*(dx**2)*(z(element)*(2.d0-xi)+z(element+1)*(1.d0+xi))
     &  /6.d0)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

      subroutine c2_quartic_coefs(nelements, x,y,yp,ypp, z)
      integer nelements
      double precision x(0:nelements),y(0:nelements),yp(0:nelements),
     &  ypp(0:0)
      double precision z(nelements)
      integer i
      double precision dx,dxnew,dxratio
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "in c2_quartic_coefs"
c     print *, "nelements = ",nelements
c     do i=0,nelements
c       print *, "x[",i,"] = ",x(i)
c     enddo
c     do i=0,nelements
c       print *, "y[",i,"] = ",y(i)
c     enddo
c     do i=0,nelements
c       print *, "yp[",i,"] = ",yp(i)
c     enddo
c     print *, "ypp[0] = ",ypp(0)
      dx=x(1)-x(0)
      z(1)=(22.d0*y(0)+10.d0*y(1)+dx*(8.d0*yp(0)-2.d0*yp(1)+dx*ypp(0)))
     &    * 0.03125d0
      do i=2,nelements
        dxnew=x(i)-x(i-1)
        dxratio=dxnew/dx
        z(i)=(-10.d0*y(i-2)*dxratio**2
     &       -2.d0*yp(i-2)*dxnew*dxratio
     &       +22.d0*y(i-1)*(1.d0-dxratio)*(1.d0+dxratio)
     &       +8.d0*yp(i-1)*dxnew*(1.+dxratio)
     &       +10.d0*y(i)
     &       -2.d0*yp(i)*dxnew)*0.03125d0
     &       +z(i-1)*dxratio**2
        dx=dxnew
      enddo
c     do i=1,nelements
c       print *, "z[",i,"] = ",z(i)
c     enddo
c     print *, "leaving c2_quartic_coefs"
c     call flush(6)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

      function c2_quartic_eval(element,nelements,pt, x,y,yp,z)
      double precision c2_quartic_eval
      integer element,nelements
      double precision pt
      double precision x(0:nelements),y(0:nelements),yp(0:nelements),
     &  z(nelements)

      double precision dx,onemxi,xi
c     integer i
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "in c2_quartic_eval"
c     print *, "element,pt",element,pt
c     do i=element,element+1
c       print *, "x[",i,"] = ",x(i)
c     enddo
c     do i=element,element+1
c       print *, "y[",i,"] = ",y(i)
c     enddo
c     do i=element,element+1
c       print *, "yp[",i,"] = ",yp(i)
c     enddo
c     print *, "z[",element+1,"] = ",z(element+1)

      dx=x(element+1)-x(element)
      xi=(pt-x(element))/dx
c     print *, "xi = ",xi
      if (xi.lt.0.d0 .or. xi.gt.1.d0) then
        print *, "in c2_quartic_eval"
        print *, "pt = ",pt
        print *, "element = ",element
        print *, "x = ",x(element),x(element+1)
        call abort()
      endif
      onemxi=1.d0-xi
      c2_quartic_eval=(1.d0-2.d0*xi)
     &  *((y(element)*(1.d0+4.d0*xi)+yp(element)*dx*xi)*onemxi**2
     &   -(y(element+1)*(5.d0-4.d0*xi)-yp(element+1)*dx*onemxi)*xi**2)
     &  +z(element+1)*16.d0*(xi*onemxi)**2
c     print *, "leaving c2_quartic_eval"
c     call flush(6)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

      function c2_spline_eval(element,nelements,order,pt,sigma,tau,
     &  x,xi,y,yp,ypp)
      double precision c2_spline_eval
      integer element,nelements,order
      double precision pt,sigma,tau
      double precision x(0:nelements),xi(0:order-4),
     &  y(0:nelements*(order-4)),yp(0:nelements),ypp(0:nelements)

      double precision dx,onemxii,prod,prod0,prod1,term1,term2,xii
      integer i,j,k
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      j=element*(order-4)
c     print *, "in c2_spline_eval"
c     print *, "element,pt",element,pt
c     do i=element,element+1
c       print *, "x[",i,"] = ",x(i)
c     enddo
c     do i=j,j+order-4
c       print *, "y[",i,"] = ",y(i)
c     enddo
c     do i=element,element+1
c       print *, "yp[",i,"] = ",yp(i)
c     enddo
c     do i=element,element+1
c       print *, "ypp[",i,"] = ",ypp(i)
c     enddo
      dx=x(element+1)-x(element)
      xii=(pt-x(element))/dx
c     print *, "xi = ",xi
      if (xii.lt.0.d0 .or. xii.gt.1.d0) then
        print *, "in c2_spline_eval"
        print *, "pt = ",pt
        print *, "element = ",element
        print *, "x = ",x(element),x(element+1)
        call abort()
      endif
      onemxii=1.d0-xii
      prod0=1.d0
      prod1=1.d0
      do i=1,order-5
        prod0=prod0*(1.d0-xii/xi(i))
        prod1=prod1*(1.d0-onemxii/xi(i))
      enddo
      term1=3.+sigma
      term2=6.+sigma*term1-0.5*tau
      c2_spline_eval=
     &  (y(j)*(1.d0+xii*(term1+xii*term2))
     &   +dx*xii*(yp(element)*(1.+xii*term1)
     &           +ypp(element)*dx*0.5*xii))*prod0*onemxii**3
     &  +(y(j+order-4)*(1.d0+onemxii*(term1+onemxii*term2))
     &   +dx*onemxii*(-yp(element+1)*(1.+onemxii*term1)
     &               +ypp(element+1)*dx*0.5*onemxii))*prod1*xii**3
      do i=1,order-5
        prod=((xii/xi(i))*(onemxii/(1.-xi(i))))**3
        do k=1,order-5
          if (k.ne.i) then
            prod=prod*(xii-xi(k))/(xi(i)-xi(k))
          endif
        enddo
        c2_spline_eval=c2_spline_eval+y(j+i)*prod
      enddo
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end
