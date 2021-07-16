c     given x and order,
c     compute Legendre polynomial f(x), and derivatives f' and f''
      subroutine legendre_polys(order,x,f,dfdx,d2fdx2)
      integer order
      double precision x
      double precision f,dfdx,d2fdx2

      double precision dfm1dx,dfp1dx,d2fm1dx2,d2fp1dx2,factor,factorm1,
     &  fm1,fp1
      integer i
c     print *, "entering legendre_polys"
c     print *, "order = ",order
c     print *, "x = ",x

      if (order.le.0) then
        f=1.d0
        dfdx=0.d0
        d2fdx2=0.d0
        return
      endif
      fm1=1.d0
      f=x
      dfm1dx=0.d0
      dfdx=1.d0
      d2fm1dx2=0.d0
      d2fdx2=0.d0
c     print *, " "
c     print *, "i = 1"
c     print *, "f,dfdx,d2fdx2 = ",f,dfdx,d2fdx2
c     print *, "fm1,dfm1dx,d2fm1dx2 = ",fm1,dfm1dx,d2fm1dx2
      do i=2,order
        factor=dble(2*i-1)/dble(i)
        factorm1=dble(i-1)/dble(i)
        fp1=factor*x*f-factorm1*fm1
        dfp1dx=factor*(x*dfdx+f)-factorm1*dfm1dx
        d2fp1dx2=factor*(x*d2fdx2+2.d0*dfdx)-factorm1*d2fm1dx2
        fm1=f
        f=fp1
        dfm1dx=dfdx
        dfdx=dfp1dx
        d2fm1dx2=d2fdx2
        d2fdx2=d2fp1dx2
c       print *, " "
c       print *, "i = ",i
c       print *, "f,dfdx,d2fdx2 = ",f,dfdx,d2fdx2
c       print *, "fm1,dfm1dx,d2fm1dx2 = ",fm1,dfm1dx,d2fm1dx2
      enddo
c     print *, "leaving legendre_polys"
c     print *, "f,dfdx,d2fdx2 = ",f,dfdx,d2fdx2
      return
      end

c***********************************************************************

c     given maximum number of iterations maxit and Legendre poly order,
c     find x between xlo and xhi
c     so that Legendre polynomial f of given order satisfies f'(x)=0
      subroutine legendre_derivative_zero(maxit,order,xlo,xhi,x)
      integer maxit,order
      double precision xhi,xlo
      double precision x
      double precision a,b,dfadx,dfbdx,dfdx,d2fadx2,d2fbdx2,d2fdx2,f,fa,
     &  fb,xnew
      integer it
c     print *, "entering legendre_derivative_zero"
c     print *, "maxit = ",maxit
c     print *, "order = ",order
c     print *, "xlo,xhi = ",xlo,xhi

      if (order.le.1) then
        x=undefind
        return
      endif
      a=xlo
      b=xhi
      call legendre_polys(order,a,fa,dfadx,d2fadx2)
      call legendre_polys(order,b,fb,dfbdx,d2fbdx2)
c     print *, "a,dfadx,d2fadx2 = ",a,dfadx,d2fadx2
c     print *, "b,dfbdx,d2fbdx2 = ",b,dfbdx,d2fbdx2
      if (dfadx*dfbdx.gt.0.d0) then
        print *, "legendre_derivative_zero not necessarily bracketed"
        call flush(6)
        call abort()
      endif
      x=0.5d0*(a+b)
      do it=1,maxit
        call legendre_polys(order,x,f,dfdx,d2fdx2)
c       print *, " "
c       print *, "a,x,b = ",a,x,b
c       print *, "dfadx,dfdx,dfbdx = ",dfadx,dfdx,dfbdx
c       print *, "x,dfdx,d2fdx2 = ",x,dfdx,d2fdx2
        if (abs(dfdx).le.roundoff*abs(f)) goto 10
        if (b-a.le.roundoff) goto 10
        if (dfdx*d2fdx2.ge.0.d0) then
          if (abs(dfdx).lt.(x-a)*abs(d2fdx2)) then
c           print *, "newton going left"
            xnew=x-dfdx/d2fdx2
          else
c           print *, "Newton going left too big: bisection"
c           if (abs(d2fdx2).ne.0.d0) then
c             print *, "x-a,dfdx/d2fdx2 = ",x-a,dfdx/d2fdx2
c           endif
            xnew=0.5d0*(a+x)
          endif
          b=x
          fb=f
          dfbdx=dfdx
          d2fbdx2=d2fdx2
        else
          xnew=x-dfdx/d2fdx2
          if (xnew.ge.b) then
c           print *, "Newton going right too big: bisection"
c           if (abs(d2fdx2).ne.0.d0) then
c             print *, "b-x,-dfdx/d2fdx2 = ",b-x,-(dfdx/d2fdx2)
c           endif
            xnew=0.5d0*(x+b)
          endif
          a=x
          fa=f
          dfadx=dfdx
          d2fadx2=d2fdx2
        endif
        x=xnew
      enddo
   10 continue
c     print *, "leaving legendre_derivative_zero"
c     print *, "x = ",x
c     call flush(6)
      return
      end

c***********************************************************************

c     given max iterations maxit and number of Lobatto nodes n,
c     find array of Lobatto nodes for all rules of order at most n
c     the nodes for the rule with i nodes with 2 <= i <= n are in
c       nodes(k+j) for 0 <=j <= i
c     where k=((i-2)*(i+1))/2
      subroutine lobatto_nodes(maxit,n,nodes)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer maxit,n
      double precision nodes((n*(n+3))/2)
      integer i,j,k
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "entering lobatto_nodes"
c     print *, "maxit = ",maxit
c     print *, "n = ",n
      k=0
      do i=0,n-1
        k=k+1
        nodes(k)=-1.d0
        k=k+i+1
        nodes(k)=1.d0
      enddo
      k=1
      do i=1,n-1
        k=k+2 ! 1st node(=-1) of lobatto rule using i zeros of legendre'
c       print *, " "
c       print *, "i = ",i
c       print *, "node indices ",k," to ",k+i+1
        do j=1,i
          k=k+1
          call legendre_derivative_zero(maxit,i+1,
     &      nodes(k-i-2),nodes(k-i-1), nodes(k))
c         if (j.eq.1 .or. j.eq.i) print *, "k = ",k
        enddo
      enddo
c     print *, "leaving lobatto_nodes"
c     k=1
c     do i=0,n-1
c       print *, "nodes = ",(nodes(k+j),j=0,i+1)
c       k=k+i+2
c     enddo
c     call flush(6)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

c***********************************************************************

c     given order and Lobatto nodes for that order,
c     find Lobatto weights, hierarchical basis functions and derivatives
      subroutine canonical(order,lobatto,
     &  basis,basis_deriv,weight,
     &  legendre,legendre_deriv)
      integer order
      double precision lobatto(0:order)

      double precision basis(0:order,0:order)
      double precision basis_deriv(0:order,0:order)
      double precision weight(0:order)

      double precision legendre(0:order),legendre_deriv(0:order)

      integer i,j
      double precision jl,j2m1,rt2j2m1,w
c     integer i
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (order.lt.1) call abort()
c     print *, "in canonical"
c     do i=0,order
c       print *, "lobatto[",i,"] = ",lobatto(i)
c     enddo

      legendre(0)=1.d0
      legendre_deriv(0)=0.d0
      legendre_deriv(1)=1.d0
      w=2.d0/(dble(order*(order+1)))
      do j=0,order
        legendre(1)=lobatto(j)
        basis(0,j)=0.5d0*(1.d0-lobatto(j))
        basis(1,j)=0.5d0*(1.d0+lobatto(j))
        basis_deriv(0,j)=-0.5d0
        basis_deriv(1,j)=0.5d0
        do i=2,order
          j2m1=2.d0*dble(i)-1
          jl=j2m1*legendre(i-1)
          rt2j2m1=1.d0/sqrt(2.d0*j2m1)
          legendre(i)=(jl*lobatto(j)-dble(i-1)*legendre(i-2))/dble(i)
          legendre_deriv(i)=jl+legendre_deriv(i-2)
          basis(i,j)=(legendre(i)-legendre(i-2))*rt2j2m1
          basis_deriv(i,j)=(legendre_deriv(i)-legendre_deriv(i-2))
     &                    *rt2j2m1
        enddo
        weight(j)=w/legendre(order)**2
      enddo

c     print *, "leaving canonical"
c     do i=0,order
c       print *, "weight[",i,"] = ",weight(i)
c     enddo
c     do i=0,order
c       print *, "basis[",i,"] = ",(basis(i,j),j=0,order)
c     enddo
c     do i=0,order
c       print *, "basis_deriv[",i,"] = ",(basis_deriv(i,j),j=0,order)
c     enddo
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

c***********************************************************************

      subroutine grid(nelements,nnodes,order,
     &  element_to_node,mesh)
      integer nelements,nnodes,order
      integer element_to_node(nelements,0:order)
      double precision mesh(0:nelements)

      integer element,i,node
      double precision dx
c     integer j
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (nnodes.ne.order*nelements+1) call abort()

      dx=1.d0/dble(nelements)
      do i=0,nelements
        mesh(i)=dble(i)*dx
      enddo

      node=1
      do element=1,nelements
        element_to_node(element,0)=node
        element_to_node(element,1)=node+order
        do i=1,order-1
          element_to_node(element,i+1)=node+i
        enddo
        node=node+order
      enddo

c     print *, "leaving grid"
c     do i=0,nelements
c       print *, "mesh[",i,"] = ",mesh(i)
c     enddo
c     do i=1,nelements
c       print *, "element_to_node[",i,"] = ",
c    &    (element_to_node(i,j),j=0,order)
c     enddo
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

c***********************************************************************

      subroutine mult(nelements,nnodes,order,
     &  basis,basis_deriv,dirichlet,element_to_node,plobatto,rlobatto,x,
     &  Ax)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer nelements,nnodes,order
      integer element_to_node(nelements,0:order)
      logical*1 dirichlet(nnodes)
      double precision 
     &  basis(0:order,0:order),
     &  basis_deriv(0:order,0:order),
     &  plobatto(nelements,0:order),
     &  rlobatto(nelements,0:order),
     &  x(nnodes)

      double precision Ax(nnodes)

      double precision bp,br
      integer element,i,j,nlobatto,node_i,node_j
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     do i=0,order
c       print *, "basis[",i,"] = ",(basis(i,j),j=0,order)
c     enddo
c     do i=0,order
c       print *, "basis_deriv[",i,"] = ",(basis_deriv(i,j),j=0,order)
c     enddo
c     do i=1,nnodes
c       print *, "dirichlet[",i,"] = ",dirichlet(i)
c     enddo
c     do i=1,nelements
c       print *, "element_to_node[",i,"] = ",
c    &    (element_to_node(i,j),j=0,order)
c     enddo
c     do i=1,nelements
c       print *, "plobatto[",i,"] = ",(plobatto(i,j),j=0,order)
c     enddo
c     do i=1,nelements
c       print *, "rlobatto[",i,"] = ",(rlobatto(i,j),j=0,order)
c     enddo
c     do i=1,nnodes
c       print *, "x[",i,"] = ",x(i)
c     enddo

      do node_i=1,nnodes
        Ax(node_i)=0.d0
      enddo
c     do i=1,nnodes
c       print *,"x,Ax[",i,"] = ",x(i),Ax(i)
c     enddo
      do element=1,nelements
c       print *, " "
c       print *, "element = ",element
        do nlobatto=0,order
c         print *, "nlobatto = ",nlobatto
          do i=0,order
            node_i=element_to_node(element,i)
            bp=basis_deriv(i,nlobatto)*plobatto(element,nlobatto)
            br=basis(i,nlobatto)*rlobatto(element,nlobatto)
c           print *, "node_i = ",node_i
c           print *, "bp,br = ",bp,br
            do j=0,order
              node_j=element_to_node(element,j)
c             print *, "x[",node_j,"] = ",x(node_j)
              Ax(node_i)=Ax(node_i)
     &          +(bp*basis_deriv(j,nlobatto)+br*basis(j,nlobatto))
     &          *x(node_j)
c             print *, "Ax[",node_i,"] = ",Ax(node_i)
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

      subroutine initialize(nelements,nnodes,order,
     &  basis,basis_deriv,element_to_node,lobatto,mesh,weight,
     &  dirichlet,plobatto,rlobatto,residual,soln)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      double precision roundoff,small,huge,undefind,pi
      common/machine/roundoff,small,huge,undefind,pi
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer nelements,nnodes,order
      integer element_to_node(nelements,0:order)
      double precision 
     &  basis(0:order,0:order),
     &  basis_deriv(0:order,0:order),
     &  lobatto(0:order),
     &  mesh(0:nelements),
     &  weight(0:order)

      logical*1 dirichlet(nnodes)
      double precision 
     &  plobatto(nelements,0:order),
     &  rlobatto(nelements,0:order),
     &  residual(nnodes),
     &  soln(nnodes)

      integer element,i,nlobatto,node_i
      double precision dx,fg,xg
c     integer j

c     bvp:
c       - d/dx( p du/dx ) + r u = f, 0 < x < 1
c       u(0)=left, u(1)=right
      double precision f,p,r,x
      double precision left,right
      p(x)=1.d0
      r(x)=0.d0
      f(x)=sin(pi*x)*pi**2
c     f(x)=1.d0

      left=0.
      right=0.
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "in initialize"
c     do i=0,order
c       print *, "basis[",i,"] = ",(basis(i,j),j=0,order)
c     enddo
c     do i=0,order
c       print *, "basis_deriv[",i,"] = ",(basis_deriv(i,j),j=0,order)
c     enddo
c     do i=1,nelements
c       print *, "element_to_node[",i,"] = ",
c    &    (element_to_node(i,j),j=0,order)
c     enddo
c     do i=0,order
c       print *, "lobatto[",i,"] = ",lobatto(i)
c     enddo
c     do i=0,nelements
c       print *, "mesh[",i,"] = ",mesh(i)
c     enddo
c     do i=0,order
c       print *, "weight[",i,"] = ",weight(i)
c     enddo
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (nnodes.ne.order*nelements+1) then
        call abort()
      endif

      do i=1,nnodes
        soln(i)=0.d0
c       let subroutine mult put inhomogeneities into initial residual
        dirichlet(i)=.false.
      enddo

c     dirichlet boundary conditions
      soln(1)=left
      soln(nnodes)=right

c     do i=1,nnodes
c       print *, "soln[",i,"] = ",soln(i)
c     enddo

      do element=1,nelements
        dx=0.5d0*(mesh(element)-mesh(element-1))
        do nlobatto=0,order
          xg=mesh(element-1)+(lobatto(nlobatto)+1.d0)*dx
          plobatto(element,nlobatto)=p(xg)*weight(nlobatto)/dx
          rlobatto(element,nlobatto)=r(xg)*weight(nlobatto)*dx
        enddo
      enddo

      call mult(nelements,nnodes,order,
     &  basis,basis_deriv,dirichlet,element_to_node,plobatto,rlobatto,
     &    soln,
     &  residual)
      do node_i=1,nnodes
        residual(node_i)=-residual(node_i)
      enddo

      do element=1,nelements
        dx=0.5d0*(mesh(element)-mesh(element-1))
c       print *, " "
c       print *, "element = ",element
        do nlobatto=0,order
          xg=mesh(element-1)+(lobatto(nlobatto)+1.d0)*dx
          fg=f(xg)*weight(nlobatto)*dx
c         print *, "fg(",xg,") = ",fg
          do i=0,order
            node_i=element_to_node(element,i)
            residual(node_i)=residual(node_i)+basis(i,nlobatto)*fg
c           print *, "basis(",i,",",nlobatto,") = ",basis(i,nlobatto)
c           print *, "residual(",node_i,") = ",residual(node_i)
          enddo
        enddo
      enddo

c     essential boundary conditions:
      dirichlet(1)=.true.
      if (dirichlet(1)) residual(1)=0.d0

      dirichlet(nnodes)=.true.
      if (dirichlet(nnodes)) residual(nnodes)=0.d0

c     do i=1,nnodes
c       print *, "residual[",i,"] = ",residual(i)
c     enddo

c     do j=1,nnodes
c       do i=1,nnodes
c         soln(i)=0.d0
c         dirichlet(i)=.false.
c       enddo
c       dirichlet(1)=.true.
c       dirichlet(nnodes)=.true.
c       soln(j)=1.d0
c       call mult(nelements,nnodes,order,
c    &    basis,basis_deriv,dirichlet,element_to_node,plobatto,rlobatto,
c    &      soln,
c    &    residual)
c       print *, " "
c       print *, "j = ",j
c       do i=1,nnodes
c         print *, "residual[",i,"] = ",residual(i)
c       enddo
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

      subroutine precg(max_cg_its,ndigit,nelements,nnodes,order,
     &  basis,basis_deriv,dirichlet,element_to_node,plobatto,rlobatto,
     &  b,x, 
     &  w)
      integer max_cg_its,ndigit,nelements,nnodes,order
      logical*1 dirichlet(nnodes)
      integer element_to_node(nelements,0:order)
      double precision 
     &  basis(0:order,0:order),
     &  basis_deriv(0:order,0:order),
     &  plobatto(nelements,0:order),
     &  rlobatto(nelements,0:order)
        
      double precision b(nnodes),x(nnodes)
      double precision w(nnodes,2)

      double precision dasum,ddot
      external dasum,ddot

      integer i,it
      double precision dif,r,s,size,t
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "max_cg_its,ndigit = ",max_cg_its,ndigit
c     print *, "nelements = ",nelements
c     do i=1,nnodes
c       print *,"x,b[",i,"] = ",x(i),b(i)
c     enddo

c     initialize variables
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
     &    basis,basis_deriv,dirichlet,element_to_node,plobatto,rlobatto,
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
        call stopit(dif,size,ndigit,max_cg_its)
c       print *, "dif = ",dif
        it=it+1
      enddo

c     print *, "leaving precg"
c     do i=1,nnodes
c       print *,"x[",i,"] = ",x(i)
c     enddo
c     call mult(nelements,nnodes,order,
c    &  basis,basis_deriv,dirichlet,element_to_node,plobatto,rlobatto,
c    &    x,
c    &  b)
c     do i=1,nnodes
c       print *,"Ax[",i,"] = ",b(i)
c     enddo
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

c***********************************************************************

      subroutine stopit(dif,size,ndigit,max_cg_its)
      double precision dif,size,e,t
      integer i,max_cg_its,ndigit
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
      else if (i.ge.max_cg_its) then
c       stopping criterion II
        e=e+2.d0
        dif=-dif
        i=0
      endif
      size=e
      return
      end

c***********************************************************************

      function approximation(nelements,nnodes,order, x,
     &  element_to_node,mesh,soln,
     &  legendre)
      double precision approximation
      integer nelements,nnodes,order
      double precision x
      integer element_to_node(nelements,0:order)
      double precision mesh(0:nelements)
      double precision soln(nnodes)
      double precision legendre(0:order)

      integer element,i
      double precision basis,dx,j2m1,xi
c     integer j
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     do i=1,nelements
c       print *, "element_to_node[",i,"] = ",
c    &    (element_to_node(i,j),j=0,order)
c     enddo
c     do i=0,nelements
c       print *, "mesh[",i,"] = ",mesh(i)
c     enddo
c     do i=1,nnodes
c       print *,"soln[",i,"] = ",soln(i)
c     enddo
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (order.lt.1) call abort()

      dx=1./dble(nelements)
      element=x/dx
      if (element.lt.nelements) element=element+1
      xi=2.d0*(x-mesh(element-1))/dx-1.d0

      approximation=0.5d0*((1.d0-xi)*soln(element_to_node(element,0))
     &                    +(1.d0+xi)*soln(element_to_node(element,1)))
      if (order.gt.1) then
        legendre(0)=1.d0
        legendre(1)=xi
        do i=2,order
          j2m1=2.d0*dble(i)-1
          legendre(i)=(j2m1*legendre(i-1)*xi-dble(i-1)*legendre(i-2))
     &               /dble(i)
          basis=(legendre(i)-legendre(i-2))/sqrt(2.d0*j2m1)
          approximation=approximation
     &                 +basis*soln(element_to_node(element,i))
        enddo
      endif
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
c     solution=0.5d0*x*(1.d0-x)
c     print *, "x = ",x
c     print *, "pi = ",pi
c     print *, "solution = ",solution
      return
      end
