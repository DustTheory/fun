      subroutine initialize_dirichlet(fi,la, mesh, matrix,rhs)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      double precision roundoff,small,huge,undefind,pi
      common/machine/roundoff,small,huge,undefind,pi
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fi,la
      double precision matrix(fi:la+1,-1:1),rhs(fi:la+1),mesh(fi-1:la+1)

      integer i
      double precision dx,x

c     bvp:
c       - d/dx( p du/dx ) + r u = f, 0 < x < 1
c       u(0)=left, u(1)=right
      double precision f,p,r
      double precision left,right
      p(x)=1.d0
      r(x)=0.d0
      f(x)=sin(pi*x)*pi**2

      left=0.
      right=0.
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (fi-1.ge.la+1) call abort()

      x=0.5d0*(mesh(fi)+mesh(fi-1))
      matrix(fi,-1)=-p(x)/(mesh(fi)-mesh(fi-1))
      do i=fi+1,la
        x=0.5d0*(mesh(i)+mesh(i-1))
        matrix(i,-1)=-p(x)/(mesh(i)-mesh(i-1))
        matrix(i-1,1)=matrix(i,-1)
      enddo
      x=0.5d0*(mesh(la+1)+mesh(la))
      matrix(la,1)=-p(x)/(mesh(la+1)-mesh(la))

      do i=fi,la
        dx=0.5d0*(mesh(i+1)-mesh(i-1))
        matrix(i,0)=-matrix(i,-1)-matrix(i,1)+r(mesh(i))*dx
        rhs(i)=f(mesh(i))*dx
      enddo

      rhs(fi)=rhs(fi)-matrix(fi,-1)*left
      rhs(la)=rhs(la)-matrix(la,1)*right
      matrix(la,1)=0.d0

c     add last equation so Dirichlet and Neumann systems same size
c       last equation involves only value at right boundary
c       last equation scaled like previous equation
      matrix(la+1,0)=matrix(la,0)
      rhs(la+1)=matrix(la+1,0)*right
      matrix(la+1,-1)=0.d0

c     call print_loc(rhs(0))
c     do i=fi,la+1
c       print *, "matrix[",i,"] = ",matrix(i,-1),matrix(i,0),matrix(i,1)
c     enddo
c     do i=fi,la+1
c       print *, "rhs[",i,"] = ",rhs(i)
c     enddo
c     print *, "leaving initialize_dirichlet"
c     call flush(6)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

c***********************************************************************

      subroutine initialize_neumann(fi,la, mesh, matrix,rhs)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      double precision roundoff,small,huge,undefind,pi
      common/machine/roundoff,small,huge,undefind,pi
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fi,la
      double precision matrix(fi:la+1,-1:1),rhs(fi:la+1),mesh(fi-1:la+1)

      integer i
      double precision dx,x

c     bvp:
c       - d/dx( p du/dx ) + r u = f, 0 < x < 1
c       u(0)=left, p(1) u'(1)=right
      double precision f,p,r
      double precision left,right
      p(x)=1.d0
      r(x)=0.d0
      f(x)=sin(pi*x)*pi**2

      left=0.
      right=-pi
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (fi-1.ge.la+1) call abort()

      x=0.5d0*(mesh(fi)+mesh(fi-1))
      matrix(fi,-1)=-p(x)/(mesh(fi)-mesh(fi-1))
      do i=fi+1,la+1
        x=0.5d0*(mesh(i)+mesh(i-1))
        matrix(i,-1)=-p(x)/(mesh(i)-mesh(i-1))
        matrix(i-1,1)=matrix(i,-1)
      enddo
c     matrix(la+1,1)=0.d0

      do i=fi,la
        dx=0.5d0*(mesh(i+1)-mesh(i-1))
        matrix(i,0)=-matrix(i,-1)-matrix(i,1)+r(mesh(i))*dx
        rhs(i)=f(mesh(i))*dx
      enddo
      dx=0.5d0*(mesh(la+1)-mesh(la))
      matrix(la+1,0)=-matrix(la+1,-1)+r(mesh(la+1))*dx
      rhs(la+1)=f(mesh(la+1))*dx

      rhs(fi)=rhs(fi)-matrix(fi,-1)*left
      rhs(la+1)=rhs(la+1)+right

c     call print_loc(rhs(0))
c     do i=fi,la+1
c       print *, "matrix[",i,"] = ",matrix(i,-1),matrix(i,0),matrix(i,1)
c     enddo
c     do i=fi,la+1
c       print *, "rhs[",i,"] = ",rhs(i)
c     enddo
c     print *, "leaving initialize_neumann"
c     call flush(6)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

c***********************************************************************

      subroutine analytical_solution(fi,la, mesh, solution)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      double precision roundoff,small,huge,undefind,pi
      common/machine/roundoff,small,huge,undefind,pi
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fi,la
      double precision mesh(fi-1:la+1)
      double precision solution(fi-1:la+1)

      integer i
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      do i=fi-1,la+1
        solution(i)=sin(pi*mesh(i))
      enddo

c     do i=fi-1,la+1
c       print *, "solution[",i,"] = ",solution(i)
c     enddo
c     print *, "leaving analytical_solution"
c     call flush(6)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

c***********************************************************************

      subroutine solve(fi,la, matrix,rhs_in_soln_out)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fi,la
      double precision matrix(fi:la+1,-1:1),rhs_in_soln_out(fi:la+1)

      integer i
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "entering solve"
c     do i=fi,la+1
c       print *, "matrix[",i,"] = ",matrix(i,-1),matrix(i,0),matrix(i,1)
c     enddo
c     do i=fi,la+1
c       print *, "rhs[",i,"] = ",rhs_in_soln_out(i)
c     enddo
c     call flush(6)

      do i=fi+1,la+1
        matrix(i,-1)=matrix(i,-1)/matrix(i-1,0)
        matrix(i,0)=matrix(i,0)-matrix(i,-1)*matrix(i-1,1)
        rhs_in_soln_out(i)=rhs_in_soln_out(i)
     &                     -matrix(i,-1)*rhs_in_soln_out(i-1)
      enddo
      rhs_in_soln_out(la+1)=rhs_in_soln_out(la+1)/matrix(la+1,0)
      do i=la,fi,-1
        rhs_in_soln_out(i)=(rhs_in_soln_out(i)
     &    -matrix(i,1)*rhs_in_soln_out(i+1))/matrix(i,0)
      enddo

c     print *, "leaving solve"
c     do i=fi,la+1
c       print *, "soln[",i,"] = ",rhs_in_soln_out(i)
c     enddo
c     call flush(6)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end
