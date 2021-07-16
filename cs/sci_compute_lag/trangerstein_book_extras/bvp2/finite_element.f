      subroutine canonical(continuity,order,
     &  basis,basis_deriv,gauss,weight)
      integer continuity,order
      double precision gauss(order)
      double precision basis(order+1,order)
      double precision basis_deriv(order+1,order)
      double precision weight(order)

      integer j
      double precision g2,g4,omg,omg2,om2g,om3g,tm3g
c     integer i
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (order.eq.1) then
c       continuous piecewise linear
        gauss(1)=0.5d0
        weight(1)=1.d0

        basis(1,1)=1.d0-gauss(1)
        basis(2,1)=gauss(1)

        basis_deriv(1,1)=-1.d0
        basis_deriv(2,1)=1.d0
      else if (order.eq.2) then
c       continuous piecewise quadratic
        gauss(1)=0.5d0*(1.d0-sqrt(1.d0/3.d0))
        gauss(2)=0.5d0*(1.d0+sqrt(1.d0/3.d0))
        weight(1)=0.5d0
        weight(2)=0.5d0

        do j=1,order
          omg=1.d0-gauss(j)
          om2g=1.d0-2.d0*gauss(j)
          g4=4.d0*gauss(j)

          basis(1,j)=om2g*omg
          basis(2,j)=g4*omg
          basis(3,j)=-(gauss(j)*om2g)

          basis_deriv(1,j)=g4-3.d0
          basis_deriv(2,j)=4.d0*om2g
          basis_deriv(3,j)=g4-1.d0
        enddo
      else if (order.eq.3) then
        gauss(1)=0.5d0*(1.d0-sqrt(0.6d0))
        gauss(2)=0.5d0
        gauss(3)=0.5d0*(1.d0+sqrt(0.6d0))
        weight(1)=5.d0/18.d0
        weight(2)=4.d0/9.d0
        weight(3)=weight(1)
        if (continuity.eq.0) then
c         continuous piecewise cubic
          do j=1,order
            omg=1.d0-gauss(j)
            om3g=1.d0-(3.d0*gauss(j))
            tm3g=1.d0+om3g

            basis(1,j)=0.5d0*om3g*tm3g*omg
            basis(2,j)=4.5d0*gauss(j)*tm3g*omg
            basis(3,j)=-(4.5d0*gauss(j)*om3g*omg)
            basis(4,j)=0.5d0*gauss(j)*om3g*tm3g

            basis_deriv(1,j)=-5.5d0+gauss(j)*(18.d0-gauss(j)*13.5d0)
            basis_deriv(2,j)=9.d0+gauss(j)*(-45.d0+gauss(j)*40.5d0)
            basis_deriv(3,j)=-4.5d0+gauss(j)*(36.d0-gauss(j)*40.5d0)
            basis_deriv(4,j)=1.d0+gauss(j)*(-9.d0+gauss(j)*13.5d0)
          enddo
        else if (continuity.eq.1) then
c         Hermite piecewise cubic
          do j=1,order
            omg=1.d0-gauss(j)
            omg2=omg**2
            g2=gauss(j)**2

c           interpolate value at 0, slope at 0, value at 1, slope at 1
            basis(1,j)=(1.d0+2.d0*gauss(j))*omg2
            basis(2,j)=gauss(j)*omg2
            basis(3,j)=(3.d0-2.d0*gauss(j))*g2
            basis(4,j)=-(omg*g2)

            basis_deriv(1,j)=-(6.d0*gauss(j)*omg)
            basis_deriv(2,j)=omg*(1.d0-3.d0*gauss(j))
            basis_deriv(3,j)=6.d0*gauss(j)*omg
            basis_deriv(4,j)=-(gauss(j)*(2.d0-3.d0*gauss(j)))
          enddo
        endif
c     else if (order.eq.4) then
c       g1=sqrt(525.d0-70.d0*sqrt(30.d0))/35.d0
c       g2=sqrt(525.d0+70.d0*sqrt(30.d0))/35.d0
c       gauss(1)=0.5d0*(1.d0-g2)
c       gauss(1)=0.5d0*(1.d0-g1)
c       gauss(1)=0.5d0*(1.d0+g1)
c       gauss(1)=0.5d0*(1.d0+g2)
c       weight(1)=(18.d0-sqrt(30.d0))/72.d0
c       weight(2)=(18.d0+sqrt(30.d0))/72.d0
c       weight(3)=weight(2)
c       weight(4)=weight(1)
c       if (continuity.eq.0) then
c       else
c         do j=1,order
c           basis(1,j)=(1.d0-2.d0*gauss(j))*(1.d0+4.d0*gauss(j))
c    &                *(1.d0-gauss(j))**2
c           basis(2,j)=gauss(j)*(2.d0*gauss(j)-1.d0)*(1.d0-gauss(j))**2
c           basis(3,j)=16.d0*(gauss(j)*(1.d0-gauss(j)))**2
c           basis(4,j)=(-1.d0+2.d0*gauss(j))*(5.d0-4.d0*gauss(j))
c    &                *gauss(j)**2
c           basis(5,j)=(1.d0-2.d0*gauss(j))*(1.d0-gauss(j))*gauss(j)**2
c           basis_deriv(1,j)=2.d0*gauss(j)*(1.d0-gauss(j))
c    &                      *(11.d0-16.d0*gauss(j))
c           basis_deriv(2,j)=(1.d0-gauss(j))
c    &                      *(1.d0+gauss(j)*(-7.d0+gauss(j)*8.d0))
c           basis_deriv(3,j)=32.d0*gauss(j)*(1.d0-gauss(j))
c    &                      *(1.d0-2.d0*gauss(j))
c           basis_deriv(4,j)=2.d0*gauss(j)*(1.d0-gauss(j))
c    &                      *(16.d0*gauss(j)-5.d0)
c           basis_deriv(5,j)=gauss(j)
c    &                      *(2.d0+gauss(j)*(-9.d0+gauss(j)*8.d0))
c         enddo
c       endif
      endif
c     do i=1,order
c       print *, "gauss[",i,"] = ",gauss(i)
c     enddo
c     do i=1,order
c       print *, "weight[",i,"] = ",weight(i)
c     enddo
c     do i=1,order+1
c       print *, "basis[",i,"] = ",(basis(i,j),j=1,order)
c     enddo
c     do i=1,order+1
c       print *, "basis_deriv[",i,"] = ",(basis_deriv(i,j),j=1,order)
c     enddo
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

c***********************************************************************

      function approximation(continuity,nelements,nnodes,order, x,
     &  element_to_node,mesh,soln)
      double precision approximation
      integer continuity,nelements,nnodes,order
      double precision x
      integer element_to_node(nelements,order+1)
      double precision mesh(nelements+1)
      double precision soln(nnodes)

      integer element
      double precision dx,omg,om3g,tm3g,xi
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      dx=1./dble(nelements)
      element=x/dx
      if (element.lt.nelements) element=element+1
      xi=(x-mesh(element))/dx
      if (order.eq.1) then
c       print *, "x = ",x,mesh(element),dx
c       print *, "element = ",element
c       print *, "xi = ",xi
c       print *, "nodes = ",element_to_node(element,1),
c    &    element_to_node(element,2)
        approximation=(1.d0-xi)*soln(element_to_node(element,1))
     &                +xi*soln(element_to_node(element,2))
      else if (order.eq.2) then
        omg=1.d0-xi
        approximation=
     &    (1.d0-2.d0*xi)*(omg*soln(element_to_node(element,1))
     &                   -xi*soln(element_to_node(element,3)))
     &   +4.d0*xi*omg*soln(element_to_node(element,2))
      else if (order.eq.3) then
        if (continuity.eq.0) then
          omg=1.d0-xi
          om3g=1.d0-(3.d0*xi)
          tm3g=1.d0+om3g
          approximation=
     &      0.5d0*om3g*tm3g*(omg*soln(element_to_node(element,1))
     &                      +xi*soln(element_to_node(element,4)))
     &     +4.5d0*xi*omg*(tm3g*soln(element_to_node(element,2))
     &                   -om3g*soln(element_to_node(element,3)))
        else
          omg=1.d0-xi
          approximation=
     &      ((1.d0+2.d0*xi)*soln(element_to_node(element,1))
     &      +xi*soln(element_to_node(element,2)))*omg**2
     &     +((3.d0-2.d0*xi)*soln(element_to_node(element,3))
     &      -omg*soln(element_to_node(element,4)))*xi**2
        endif
      endif
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

c***********************************************************************

      subroutine grid(continuity,nelements,nnodes,order,
     &  element_to_node,mesh)
      integer continuity,nelements,nnodes,order
      integer element_to_node(nelements,order+1)
      double precision mesh(nelements+1)

      integer element,e2,e3,i
      double precision dx
c     integer j
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (nnodes.ne.(order-continuity)*nelements+continuity+1) then
        call abort()
      endif

      dx=1.d0/dble(nelements)
      do i=1,nelements+1
        mesh(i)=dble(i-1)*dx
      enddo

      if (continuity.ge.order) call abort()
      if (order.eq.1) then
        if (continuity.ne.0) call abort()
        do element=1,nelements
          element_to_node(element,1)=element
          element_to_node(element,2)=element+1
        enddo
      else if (order.eq.2) then
        if (continuity.ne.0) call abort()
        do element=1,nelements
          e2=2*element
          element_to_node(element,1)=e2-1
          element_to_node(element,2)=e2
          element_to_node(element,3)=e2+1
        enddo
        e2=2*nelements
      else if (order.eq.3) then
        if (continuity.eq.0) then
          do element=1,nelements
            e3=3*element
            element_to_node(element,1)=e3-2
            element_to_node(element,2)=e3-1
            element_to_node(element,3)=e3
            element_to_node(element,4)=e3+1
          enddo
          e3=3*nelements
        else
          if (continuity.ne.1) call abort()
          do element=1,nelements
            e2=2*element
            element_to_node(element,1)=e2-1
            element_to_node(element,2)=e2
            element_to_node(element,3)=e2+1
            element_to_node(element,4)=e2+2
          enddo
          e2=2*nelements
        endif
      endif

c     do i=1,nelements+1
c       print *, "mesh[",i,"] = ",mesh(i)
c     enddo
c     do i=1,nelements
c       print *, "element_to_node[",i,"] = ",
c    &    (element_to_node(i,j),j=1,order+1)
c     enddo
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

c***********************************************************************

      function solution(x)
      double precision solution
      double precision x
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      double precision roundoff,small,huge,undefind,pi
      common/machine/roundoff,small,huge,undefind,pi
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      solution=sin(pi*x)
c     print *, "x = ",x
c     print *, "pi = ",pi
c     print *, "solution = ",solution
      return
      end

c***********************************************************************

      subroutine mult(nelements,nnodes,order,
     &  basis,basis_deriv,dirichlet,element_to_node,pgauss,rgauss,x,
     &  Ax)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer nelements,nnodes,order
      integer element_to_node(nelements,order+1)
      logical*1 dirichlet(nnodes)
      double precision basis(order+1,order),basis_deriv(order+1,order),
     &  pgauss(nelements,order),rgauss(nelements,order),x(nnodes)

      double precision Ax(nnodes)

      double precision bp,br
      integer element,i,j,ngauss,node_i,node_j
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     do i=1,order+1
c       print *, "basis[",i,"] = ",(basis(i,j),j=1,order)
c     enddo
c     do i=1,order+1
c       print *, "basis_deriv[",i,"] = ",(basis_deriv(i,j),j=1,order)
c     enddo
c     do i=1,nnodes
c       print *, "dirichlet[",i,"] = ",dirichlet(i)
c     enddo
c     do i=1,nelements
c       print *, "element_to_node[",i,"] = ",
c    &    (element_to_node(i,j),j=1,order+1)
c     enddo
c     do i=1,nelements
c       print *, "pgauss[",i,"] = ",(pgauss(i,j),j=1,order)
c     enddo
c     do i=1,nelements
c       print *, "rgauss[",i,"] = ",(rgauss(i,j),j=1,order)
c     enddo
c     do i=1,nnodes
c       print *, "x[",i,"] = ",x(i)
c     enddo

      do node_i=1,nnodes
        Ax(node_i)=0.d0
      enddo
      do element=1,nelements
        do ngauss=1,order
          do i=1,order+1
            node_i=element_to_node(element,i)
            bp=basis_deriv(i,ngauss)*pgauss(element,ngauss)
            br=basis(i,ngauss)*rgauss(element,ngauss)
            do j=1,order+1
              node_j=element_to_node(element,j)
              Ax(node_i)=Ax(node_i)
     &          +(bp*basis_deriv(j,ngauss)+br*basis(j,ngauss))*x(node_j)
            enddo
          enddo
        enddo
      enddo
      do node_i=1,nnodes
        if (dirichlet(node_i)) Ax(node_i)=0.d0
      enddo
c     do i=1,nnodes
c       print *,"x,Ax[",i,"] = ",x(i),Ax(i)
c     enddo
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

c***********************************************************************

      subroutine initialize(continuity,nelements,nnodes,order,
     &  basis,basis_deriv,element_to_node,gauss,mesh,weight,
     &  dirichlet,pgauss,rgauss,residual,soln)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      double precision roundoff,small,huge,undefind,pi
      common/machine/roundoff,small,huge,undefind,pi
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer continuity,nelements,nnodes,order
      integer element_to_node(nelements,order+1)
      double precision basis(order+1,order),basis_deriv(order+1,order),
     &  gauss(order),mesh(nelements+1),weight(order)

      logical*1 dirichlet(nnodes)
      double precision pgauss(nelements,order),rgauss(nelements,order),
     &  residual(nnodes),soln(nnodes)

      integer element,i,j,ngauss,node_i
      double precision dx,fg,xg

      double precision solution
      external solution

c     bvp:
c       - d/dx( p du/dx ) + r u = f, 0 < x < 1
c       u(0)=left, u(1)=right
      double precision f,p,r,x
      double precision left,right
      p(x)=1.d0
      r(x)=0.d0
      f(x)=sin(pi*x)*pi**2

      left=0.
      right=0.
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (nnodes.ne.(order-continuity)*nelements+continuity+1) then
        call abort()
      endif

      do i=1,nnodes
        soln(i)=0.d0
c       let subroutine mult put inhomogeneities into initial residual
        dirichlet(i)=.false.
      enddo

c     dirichlet boundary conditions
      j=1
      soln(j)=left
      j=nnodes-continuity
      soln(j)=right

c     do i=1,nnodes
c       print *, "soln[",i,"] = ",soln(i)
c     enddo

      do element=1,nelements
        dx=mesh(element+1)-mesh(element)
        do ngauss=1,order
          xg=mesh(element)+gauss(ngauss)*dx
          pgauss(element,ngauss)=p(xg)*weight(ngauss)/dx
          rgauss(element,ngauss)=r(xg)*weight(ngauss)*dx
        enddo
      enddo

      call mult(nelements,nnodes,order,
     &  basis,basis_deriv,dirichlet,element_to_node,pgauss,rgauss,soln,
     &  residual)
      do node_i=1,nnodes
        residual(node_i)=-residual(node_i)
      enddo

      do element=1,nelements
        dx=mesh(element+1)-mesh(element)
        do ngauss=1,order
          xg=mesh(element)+gauss(ngauss)*dx
          fg=f(xg)*weight(ngauss)*dx
          do i=1,order+1
            node_i=element_to_node(element,i)
            residual(node_i)=residual(node_i)+basis(i,ngauss)*fg
          enddo
        enddo
      enddo

c     essential boundary conditions:
      j=1
      dirichlet(j)=.true.
      if (dirichlet(j)) residual(j)=0.d0

      j=nnodes-continuity
      dirichlet(j)=.true.
      if (dirichlet(j)) residual(j)=0.d0

c     do i=1,nnodes
c       print *, "residual[",i,"] = ",residual(i)
c     enddo
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

c***********************************************************************

      subroutine pre(nnodes,x, Px)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer nnodes
      double precision x(nnodes)
      double precision Px(nnodes)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call dcopy(nnodes,x,1,Px,1)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

c***********************************************************************

      subroutine precg(limit,ndigit,nelements,nnodes,order,
     &  basis,basis_deriv,dirichlet,element_to_node,pgauss,rgauss,
     &  b,x, 
     &  w)
      integer limit,ndigit,nelements,nnodes,order
      logical*1 dirichlet(nnodes)
      integer element_to_node(nelements,order+1)
      double precision basis(order+1,order),basis_deriv(order+1,order),
     &  pgauss(nelements,order),rgauss(nelements,order)
        
      double precision b(nnodes),x(nnodes)
      double precision w(nnodes,2)

      double precision dasum,ddot
      external dasum,ddot

      integer i,it
      double precision dif,r,s,size,t
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     initialize variables
c     print *, "limit,ndigit = ",limit,ndigit
c     print *, "nelements = ",nelements
c     do i=1,nnodes
c       print *,"x,b[",i,"] = ",x(i),b(i)
c     enddo

c     call mult(nnodes,order,
c    &  dirichlet,matrix,x, 
c    &  w)
c     do i=1,nnodes
c       print *,"Ax[",i,"] = ",w(i,1)
c     enddo

c     assume that b contains the initial residual
      do i = 1,nnodes
c       w(i,1) = b(i) - w(i,1)
        w(i,1) = b(i)
      enddo
c     do i=1,nnodes
c       print *,"b-Ax[",i,"] = ",w(i,1)
c     enddo

      call pre(nnodes,w, w(1,2))
c     do i=1,nnodes
c       print *,"P*(b-Ax)[",i,"] = ",w(i,2)
c     enddo

      s = ddot(nnodes,w(1,1),1,w(1,2),1)
c     print *, "r . P r = ",s

      if ( s .le. 0.d0 ) return
c     start iterations 
      dif=1.d0
      it=0
      do while (dif.gt.0.d0)
        call mult(nelements,nnodes,order,
     &    basis,basis_deriv,dirichlet,element_to_node,pgauss,rgauss,
     &      w(1,2), 
     &    b)
c       print *, " "
c       print *, "it = ",it
c       do i=1,nnodes
c         print *,"Ap[",i,"] = ",b(i)
c       enddo

        t = ddot(nnodes,w(1,2),1,b,1)
c       print *, "p . A p = ",t

        if ( t .le. 0.d0 ) return
        r = s/t
        dif = 0.d0
        size = 0.d0
        call daxpy(nnodes,r,w(1,2),1,x,1)
c       do i=1,nnodes
c         print *,"x[",i,"] = ",x(i)
c       enddo

        call daxpy(nnodes,-r,b,1,w(1,1),1)
c       do i=1,nnodes
c         print *,"r[",i,"] = ",w(i,1)
c       enddo

        dif=r*dasum(nnodes,w(1,2),1)
        size=dasum(nnodes,x,1)
c       print *, "dif = ",dif
c       print *, "size = ",size

        call pre(nnodes,w, b)
c       do i=1,nnodes
c         print *,"P r[",i,"] = ",b(i)
c       enddo

        t = ddot(nnodes,b,1,w(1,1),1)
c       print *, "r . P r = ",t
       
        r = t/s
        do i = 1,nnodes
          w(i,2) = b(i) + r*w(i,2)
        enddo
c       do i=1,nnodes
c         print *,"p[",i,"] = ",w(i,2)
c       enddo

        s = t
        call stopit(dif,size,ndigit,limit)
c       print *, "dif = ",dif
        it=it+1
      enddo
c     do i=1,nnodes
c       print *,"x[",i,"] = ",x(i)
c     enddo
c     call mult(nelements,nnodes,order,
c    &  basis,basis_deriv,dirichlet,element_to_node,pgauss,rgauss,
c    &    x,
c    &  b)
c     do i=1,nnodes
c       print *,"Ax[",i,"] = ",b(i)
c     enddo
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

c***********************************************************************

      subroutine stopit(dif,size,ndigit,limit)
      double precision dif,size,e,t
      integer i,limit,ndigit
      data i/0/
      save i,t
      dif = abs(dif)
      size = abs(size)
c     initialization during first iteration 
      if ( i .le. 0 ) t = 10.**(-ndigit)
      i = i + 1
      e = 3.d0*dble(i)
c     print *, "i,t,e = ",i,t,e

      if (dif.le.t*size) then
c       stopping criterion I 
        e=e+1.d0
        dif=-dif
        i=0
      else if (i.ge.limit) then
c       stopping criterion II
        e=e+2.d0
        dif=-dif
        i=0
      endif
      size=e
      return
      end
