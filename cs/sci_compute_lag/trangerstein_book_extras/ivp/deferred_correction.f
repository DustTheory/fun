c     spectral deferred correction needs to compute integrals of
c       Lagrange interpolation polynomials
c     given Lobatto rule with order+2 nodes, arrays of nodes and weights
c     compute
c       integrals(i,j) = int_(-1)^t_i Lagrange polynomial_j(t) d t
      subroutine lagrange_polynomial_integrals(order,nodes,weights,
     &  integrals)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer order
      double precision nodes(0:order),weights(0:order)
      double precision integrals(1:order,0:order)
      integer i,j,k,n
      double precision denom,dt,num,sum,t,tbar
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "entering lagrange_polynomial_integrals"
c     print *, "order = ",order
c     print *, "nodes = ",(nodes(i),i=0,order)
      do k=1,order ! integrate lagrange poly from t_0 to t_k
        dt=0.5d0*(nodes(k)-nodes(0))
        tbar=0.5d0*(nodes(k)+nodes(0))
        do j=0,order ! lagrange poly omits node j
          sum=0.d0
          denom=1.d0
          do i=0,order
            if (i.ne.j) denom=denom*(nodes(j)-nodes(i))
          enddo
          do n=0,order ! quadrature point in lobatto rule
            t=nodes(n)*dt+tbar ! mapped quadrature point in (t_0,t_k)
            num=1.d0
            do i=0,order
              if (i.ne.j) num=num*(t-nodes(i))
            enddo
            sum=sum+weights(n)*num/denom
          enddo
          integrals(k,j)=sum*dt
        enddo
      enddo
c     print *, "leaving lagrange_polynomial_integrals"
c     do k=1,order
c       print *, " integrals from",nodes(0)," to ",nodes(k)," = ",
c    &    (integrals(k,j),j=0,order)
c     enddo
c     call flush(6)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

c***********************************************************************

      subroutine deferred_correction(order,nodes,integrals,dt,t,y,f,
     &  farray,fnewarray,tarray,yarray,ynewarray)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer order
      double precision nodes(0:order),integrals(1:order,0:order),
     &  dt,t,y,f
      external f
      double precision farray(0:order),fnewarray(0:order),
     &  tarray(0:order),yarray(0:order),ynewarray(0:order)
      integer i
      double precision d,eold,enew
      integer j,k
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "entering deferred_correction"
c     print *, "order = ",order
c     print *, "t,dt = ",t,dt
c     print *, "y = ",y
c     print *, "nodes = ",(nodes(i),i=0,order)
c     do k=1,order
c       print *, " integrals from",nodes(0)," to ",nodes(k)," = ",
c    &    (integrals(k,j),j=0,order)
c     enddo
c     call flush(6)
      yarray(0)=y
      tarray(0)=t
      farray(0)=f(tarray(0),yarray(0))
      fnewarray(0)=farray(0)
      ynewarray(0)=yarray(0)
      do i=0,order-1
        tarray(i+1)=t+0.5d0*dt*(nodes(i+1)+1.d0)
        yarray(i+1)=yarray(i)+(tarray(i+1)-tarray(i))*farray(i)
        farray(i+1)=f(tarray(i+1),yarray(i+1))
c       print *, "forward euler y(",tarray(i+1),") = ",yarray(i+1)
      enddo
      do j=1,order-1
        eold=0.d0
        d=0.d0
c       print *, "correction ",j
        do i=0,order-1
          enew=0.d0
          do k=0,order
            enew=enew+farray(k)*integrals(i+1,k)
          enddo
          enew=enew*0.5d0*dt+(yarray(0)-yarray(i+1))
          d=d+(tarray(i+1)-tarray(i))*(fnewarray(i)-farray(i))
     &     +(enew-eold)
          eold=enew
          ynewarray(i+1)=yarray(i+1)+d
          fnewarray(i+1)=f(tarray(i+1),ynewarray(i+1))
c         print *, "at t = ",tarray(i+1)," d,e,y,f = ",d,enew,
c    &      ynewarray(i+1),fnewarray(i+1)
        enddo
        do i=1,order
          yarray(i)=ynewarray(i)
          farray(i)=fnewarray(i)
        enddo
c       print *, "at t = ",tarray(order)," y = ",yarray(order)
      enddo
      y=yarray(order)
c     print *, "leaving deferred_correction"
c     print *, "y = ",y
c     call flush(6)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

      subroutine deferred_correctione(order,nodes,integrals,dt,t,y,f,
     &  farray,fnewarray,tarray,yarray,ynewarray)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer order
      double precision nodes(0:order),integrals(1:order,0:order),
     &  dt,t
      double complex y,f
      external f
      double precision tarray(0:order)
      double complex farray(0:order),fnewarray(0:order),yarray(0:order),
     &  ynewarray(0:order)
      integer i
      double complex d,eold,enew
      integer j,k
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "entering deferred_correctionc"
c     print *, "order = ",order
c     print *, "t,dt = ",t,dt
c     print *, "y = ",y
c     print *, "nodes = ",(nodes(i),i=0,order)
c     do k=1,order
c       print *, " integrals from",nodes(0)," to ",nodes(k)," = ",
c    &    (integrals(k,j),j=0,order)
c     enddo
c     call flush(6)
      yarray(0)=y
      tarray(0)=t
      farray(0)=f(tarray(0),yarray(0))
      fnewarray(0)=farray(0)
      ynewarray(0)=yarray(0)
      do i=0,order-1
        tarray(i+1)=t+0.5d0*dt*(nodes(i+1)+1.d0)
        yarray(i+1)=yarray(i)+(tarray(i+1)-tarray(i))*farray(i)
        farray(i+1)=f(tarray(i+1),yarray(i+1))
c       print *, "forward euler y(",tarray(i+1),") = ",yarray(i+1)
      enddo
      do j=1,order-1
        eold=cmplx(0.d0,0.d0)
        d=cmplx(0.d0,0.d0)
c       print *, "correction ",j
        do i=0,order-1
          enew=cmplx(0.d0,0.d0)
          do k=0,order
            enew=enew+farray(k)*integrals(i+1,k)
          enddo
          enew=enew*0.5d0*dt+(yarray(0)-yarray(i+1))
          d=d+(tarray(i+1)-tarray(i))*(fnewarray(i)-farray(i))
     &     +(enew-eold)
          eold=enew
          ynewarray(i+1)=yarray(i+1)+d
          fnewarray(i+1)=f(tarray(i+1),ynewarray(i+1))
c         print *, "at t = ",tarray(i+1)," d,e,y,f = ",d,enew,
c    &      ynewarray(i+1),fnewarray(i+1)
        enddo
        do i=1,order
          yarray(i)=ynewarray(i)
          farray(i)=fnewarray(i)
        enddo
c       print *, "at t = ",tarray(order)," y = ",yarray(order)
      enddo
      y=yarray(order)
c     print *, "leaving deferred_correction"
c     print *, "y = ",y
c     call flush(6)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

      subroutine deferred_correctionbe(order,nodes,integrals,dt,t,y,f,
     &  dfdy, farray,fnewarray,tarray,yarray,ynewarray)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer order
      double precision nodes(0:order),integrals(1:order,0:order),
     &  dt,t
      double complex y,f,dfdy
      external dfdy,f
      double precision tarray(0:order)
      double complex farray(0:order),fnewarray(0:order),yarray(0:order),
     &  ynewarray(0:order)
      double complex d,dfvaldy,dnew,eold,enew,fval
      double precision tdif
      integer i,it,j,k
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "entering deferred_correctionbe"
c     print *, "order = ",order
c     print *, "t,dt = ",t,dt
c     print *, "y = ",y
c     print *, "nodes = ",(nodes(i),i=0,order)
c     do k=1,order
c       print *, " integrals from",nodes(0)," to ",nodes(k)," = ",
c    &    (integrals(k,j),j=0,order)
c     enddo
c     call flush(6)
      yarray(0)=y
      tarray(0)=t
      farray(0)=f(tarray(0),yarray(0))
      fnewarray(0)=farray(0)
      ynewarray(0)=yarray(0)
      do i=0,order-1
        tarray(i+1)=t+0.5d0*dt*(nodes(i+1)+1.d0)
        yarray(i+1)=yarray(i)
        tdif=tarray(i+1)-tarray(i)
        do it=0,5
          farray(i+1)=f(tarray(i+1),yarray(i+1))
          fval=yarray(i+1)-yarray(i)-tdif*farray(i+1)
          dfvaldy=1.-tdif*dfdy(tarray(i+1),yarray(i+1))
          yarray(i+1)=yarray(i+1)-fval/dfvaldy
        enddo
        farray(i+1)=f(tarray(i+1),yarray(i+1))
c       print *, "forward euler y(",tarray(i+1),") = ",yarray(i+1)
      enddo
      do j=1,order-1
        eold=cmplx(0.d0,0.d0)
        d=cmplx(0.d0,0.d0)
c       print *, "correction ",j
        do i=0,order-1
          enew=cmplx(0.d0,0.d0)
          do k=0,order
            enew=enew+farray(k)*integrals(i+1,k)
          enddo
          enew=enew*0.5d0*dt+(yarray(0)-yarray(i+1))
          dnew=d
          tdif=tarray(i+1)-tarray(i)
          do it=0,5
            ynewarray(i+1)=yarray(i+1)+dnew
            fnewarray(i+1)=f(tarray(i+1),ynewarray(i+1))
            fval=dnew-d-tdif*(fnewarray(i+1)-farray(i+1))-(enew-eold)
            dfvaldy=1.-tdif*dfdy(tarray(i+1),ynewarray(i+1))
            dnew=dnew-fval/dfvaldy
          enddo
          eold=enew
          d=dnew
          ynewarray(i+1)=yarray(i+1)+d
          fnewarray(i+1)=f(tarray(i+1),ynewarray(i+1))
c         print *, "at t = ",tarray(i+1)," d,e,y,f = ",d,enew,
c    &      ynewarray(i+1),fnewarray(i+1)
        enddo
        do i=1,order
          yarray(i)=ynewarray(i)
          farray(i)=fnewarray(i)
        enddo
c       print *, "at t = ",tarray(order)," y = ",yarray(order)
      enddo
      y=yarray(order)
c     print *, "leaving deferred_correction"
c     print *, "y = ",y
c     call flush(6)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end
