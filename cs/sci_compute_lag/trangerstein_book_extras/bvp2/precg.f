      subroutine precg(limit,ndigit,nnodes,order,
     &  dirichlet,matrix,
     &  b,x, 
     &  w)
c input
      integer limit,ndigit,nnodes,order
      logical*1 dirichlet(nnodes)
      double precision matrix(nnodes,-order:order)
c in-out
      double precision b(nnodes),x(nnodes)
c work
      double precision w(nnodes,2)
c local
      integer i
      double precision dif,r,s,size,t
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     initialize variables
      call mult(nnodes,order,
     &  dirichlet,matrix,x, 
     &  w)
      do i = 1,nnodes
        w(i,1) = b(i) - w(i,1)
      enddo
      call pre(nnodes,w, w(1,2))
      s = ddot(nnodes,w(1,1),1,w(1,2),1)
      if ( s .le. 0.d0 ) return
c     start iterations 
      dif=1.d0
      do while (dif.gt.0.d0)
        call mult(nnodes,order,
     &    dirichlet,matrix,b, 
     &    w(1,2))
        t = ddot(nnodes,w(1,2),1,b,1)
        if ( t .le. 0.d0 ) call abort()
        r = s/t
        dif = 0.d0
        size = 0.d0
        call daxpy(nnodes,r,w(1,2),1,x,1)
        call daxpy(nnodes,-r,b,1,w(1,1),1)
        dif=r*dasum(nnodes,w(1,2),1)
        size=dasum(nnodes,x,1)
        call pre(nnodes,w, b)
        t = ddot(nnodes,b,1,w(1,1),1)
        r = t/s
        do i = 1,nnodes
          w(i,2) = b(i) + r*w(i,2)
        enddo
        s = t
        call stopit(dif,size,ndigit,limit)
      enddo
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

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
      if (dif.le.t*size) then
c       stopping criterion I 
        e=e+1.d0
      else if (i.ge.limit) then
c       stopping criterion II
        e=e+2.d0
      endif
      if (dif.le.t*size .or. i.lt.limit) then
        dif=-dif
        i=0
      endif
      size=e
      return
      end
