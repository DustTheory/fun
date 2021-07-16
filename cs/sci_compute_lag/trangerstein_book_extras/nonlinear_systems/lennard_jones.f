      program lennard_jones
      integer atoms,itmax,lwa,n
      parameter (atoms=2)
      parameter (n=1) ! max(1,3*atoms-6)
      parameter (lwa=19) ! (n*(3*n+13))/2
      parameter (itmax=1)

c     parameter (atoms=3)
c     parameter (n=3) ! max(1,3*atoms-6)
c     parameter (lwa=33) ! (n*(3*n+13))/2
c     parameter (itmax=20)
      double precision a(3,atoms),x(n),f(n),wa(lwa)
      integer i,info,it,k
      double precision low,high,len,tol
      external lj,ljgrad
      double precision lj

c     n=max(1,3*atoms-6)
c     lwa=(n*(3*n+13))/2

c     low=0.1d0
c     high=1.2d0
c     len=high-low
c     do i=1,n
c       x(i)=low+len*rand()
c     enddo
      x(1)=2.d0**(1.d0/6.d0)
      if (atoms.ge.3) then
        x(2)=0.5d0*x(1)
        x(3)=0.5d0*x(1)*sqrt(3.d0)
      endif
      if (atoms.ge.4) then
        x(4)=x(2)
        x(5)=x(1)/sqrt(12.d0)
        x(6)=x(1)*sqrt(2.d0/3.d0)
      endif
      do i=1,3
        a(i,1)=0.d0
      enddo
      a(1,2)=x(1)
      a(2,2)=0.d0
      a(3,2)=0.d0
      if (atoms.ge.3) then
        a(1,3)=x(2)
        a(2,3)=x(3)
        a(3,3)=0.d0
      endif
      do k=4,atoms
        do i=1,3
          a(i,k)=x(i+3*k-9)
        enddo
      enddo
      do k=1,atoms
        print *, "atom[",k,"] = ",(a(i,k),i=1,3)
      enddo
      print *, "optimal energy = ",lj(x)

      do it=1,itmax
        x(1)=2.d0**(1.d0/6.d0)
        if (atoms.ge.3) then
c         x(2)=0.5d0*x(1)
c         x(3)=0.5d0*x(1)*sqrt(3.d0)
          x(2)=x(1)*rand()
          x(3)=x(1)*rand()
        endif
        if (atoms.ge.4) then
          x(4)=x(2)
          x(5)=x(1)/sqrt(12.d0)
          x(6)=x(1)*sqrt(2.d0/3.d0)
        endif
        tol=1.d-12

        call hybrd1(ljgrad,n,x,f,tol,info,wa,lwa)

        do i=1,3
          a(i,1)=0.d0
        enddo
        a(1,2)=x(1)
        a(2,2)=0.d0
        a(3,2)=0.d0
        if (atoms.ge.3) then
          a(1,3)=x(2)
          a(2,3)=x(3)
          a(3,3)=0.d0
        endif
        do k=4,atoms
          do i=1,3
            a(i,k)=x(i+3*k-9)
          enddo
        enddo

        print *, "info = ",info
        print *, "energy = ",lj(x)
        do k=1,atoms
          print *, "atom[",k,"] = ",(a(i,k),i=1,3)
        enddo
      enddo

      stop
      end

      double precision function lj(x)
      integer atoms,n
      parameter (atoms=2)
      parameter (n=1) ! max(1,3*atoms-6)
c     parameter (atoms=3)
c     parameter (n=3) ! max(1,3*atoms-6)
      double precision a(3,atoms),x(n)
      integer i,j
      double precision d,t

      do i=1,3
        a(i,1)=0.d0
      enddo
      a(1,2)=x(1)
      a(2,2)=0.d0
      a(3,2)=0.d0
      if (atoms.ge.3) then
        a(1,3)=x(2)
        a(2,3)=x(3)
        a(3,3)=0.
      endif
      do j=4,atoms
        do i=1,3
          a(i,j)=x(i+3*j-9)
        enddo
      enddo

      lj=0.d0
      do i=1,atoms-1
        do j=i+1,atoms
          d=(a(1,i)-a(1,j))**2+(a(2,i)-a(2,j))**2+(a(3,i)-a(3,j))**2
          t=d**(-6)-d**(-3)
c         print *, "dist[",i,",",j,"] = ",d," t = ",t
          lj=lj+d**(-6)-d**(-3)
        enddo
      enddo
c     print *, "lj = ",lj

      return
      end

      subroutine ljgrad(n,x,fvec,iflag)
      integer n,iflag
      double precision x(n),fvec(n)
      integer atoms
      parameter (atoms=2)
c     parameter (atoms=3)
      double precision a(3,atoms)
      integer i,j,k
      double precision d,t

      do i=1,3
        a(i,1)=0.d0
      enddo
      a(1,2)=x(1)
      a(2,2)=0.d0
      a(3,2)=0.d0
      if (atoms.ge.3) then
        a(1,3)=x(2)
        a(2,3)=x(3)
        a(3,3)=0.d0
      endif
      do j=4,atoms
        do i=1,3
          a(i,j)=x(i+3*j-9)
        enddo
      enddo

      i=1
      j=2
      s=0.d0
      do k=1,atoms
        if (k .ne. j) then
          d=(a(1,j)-a(1,k))**2+(a(2,j)-a(2,k))**2+(a(3,j)-a(3,k))**2
          s=s+(a(i,j)-a(i,k))*(d**3-2.d0)/d**7
        endif
        fvec(1)=6.d0*s
      enddo

      j=3
      fvec(2)=0.d0
      fvec(3)=0.d0
      do k=1,atoms
        if (k.ne.j) then
          d=(a(1,j)-a(1,k))**2+(a(2,j)-a(2,k))**2+(a(3,j)-a(3,k))**2
          t=(d**3-2.d0)/d**7
          do i=1,2
            fvec(i+1)=fvec(i+1)+(a(i,j)-a(i,k))*t
          enddo
        endif
      enddo
      fvec(2)=6.d0*fvec(2)
      fvec(3)=6.d0*fvec(3)

      do j=4,atoms
        do i=1,3
          fvec(i+3*j-9)=0.d0
        enddo
        do k=1,atoms
          if (k.ne.j) then
            d=(a(1,j)-a(1,k))**2+(a(2,j)-a(2,k))**2+(a(3,j)-a(3,k))**2
            t=(d**3-2.d0)/d**7
            do i=1,2
              fvec(i+3*j-9)=fvec(i+3*j-9)+(a(i,j)-a(i,k))*t
            enddo
          endif
        enddo
        do i=1,3
          fvec(i+3*j-9)=6.d0*fvec(i+3*j-9)
        enddo
      enddo

      return
      end
