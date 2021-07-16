      subroutine matrix_multiply(fi,la,ifirst,ilast,
     &  a,x, 
     &  ax)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fi,la,ifirst,ilast
      double precision 
     &  a(fi:la,-1:1),
     &  x(fi:la)
      double precision 
     &  ax(fi:la)

      integer ic
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     dirichlet bc at left
      ic=ifirst
      ax(ic)=a(ic,0)*x(ic)+a(ic,1)*x(ic+1)
      do ic=ifirst+1,ilast-1
        ax(ic)=a(ic,-1)*x(ic-1)+a(ic,0)*x(ic)+a(ic,1)*x(ic+1)
      enddo
c     dirichlet bc at right
      ic=ilast
      ax(ic)=a(ic,-1)*x(ic-1)+a(ic,0)*x(ic)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

      subroutine compute_residual(fi,la,ifirst,ilast,
     &  a,b,x, 
     &  ax_b)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fi,la,ifirst,ilast
      double precision 
     &  a(fi:la,-1:1),
     &  b(fi:la),
     &  x(fi:la)
      double precision 
     &  ax_b(fi:la)

      integer ic
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ic=ifirst
      ax_b(ic)=a(ic,0)*x(ic)+a(ic,1)*x(ic+1)-b(ic)
      do ic=ifirst+1,ilast-1
        ax_b(ic)=a(ic,-1)*x(ic-1)+a(ic,0)*x(ic)+a(ic,1)*x(ic+1)-b(ic)
      enddo
      ic=ilast
      ax_b(ic)=a(ic,-1)*x(ic-1)+a(ic,0)*x(ic)-b(ic)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 


      subroutine richardson(fi,la,ifirst,ilast, mu,
     &  matrix,rhs,
     &  solution,
     &  residual)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      double precision roundoff,small,huge,undefind
      common/machine/roundoff,small,huge,undefind
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fi,la,ifirst,ilast
      double precision mu
      double precision matrix(fi:la,-1:1),rhs(fi:la)
      double precision solution(fi:la)
      double precision residual(fi:la)

      integer ic
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "in richardson, fi,la = ",fi,la
c     print *, "ifirst,ilast = ",ifirst,ilast
c     print *, "mu = ",mu
c     do ic=ifirst,ilast
c       print *, "matrix[",ic,"]= ",matrix(ic,-1),matrix(ic,0),
c    &    matrix(ic,1)
c     enddo
c     do ic=ifirst,ilast
c       print *, "rhs[",ic,"]= ",rhs(ic)
c     enddo
c     do ic=ifirst,ilast
c       print *, "solution[",ic,"]= ",solution(ic)
c     enddo

      call compute_residual(fi,la,ifirst,ilast,matrix,rhs,solution,
     &  residual)
c     do ic=ifirst,ilast
c       print *, "residual[",ic,"]= ",residual(ic)
c     enddo

      do ic=ifirst,ilast
        solution(ic)=solution(ic)-residual(ic)/mu
      enddo
c     do ic=ifirst,ilast
c       print *, "solution[",ic,"]= ",solution(ic)
c     enddo
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

      subroutine jacobi(fi,la,ifirst,ilast,
     &  matrix,rhs,
     &  solution,
     &  residual)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fi,la,ifirst,ilast
      double precision matrix(fi:la,-1:1),rhs(fi:la)
      double precision solution(fi:la)
      double precision residual(fi:la)

      integer ic
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call compute_residual(fi,la,ifirst,ilast,matrix,rhs,solution,
     &  residual)
      do ic=ifirst,ilast
        solution(ic)=solution(ic)-residual(ic)/matrix(ic,0)
      enddo
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

      subroutine jacobi_preconditioner(fi,la,ifirst,ilast,
     &  matrix,rhs,
     &  solution)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fi,la,ifirst,ilast
      double precision matrix(fi:la,-1:1),rhs(fi:la)
      double precision solution(fi:la)

      integer ic
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      do ic=ifirst,ilast
        solution(ic)=rhs(ic)/matrix(ic,0)
      enddo
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

      subroutine jacobi_omega(fi,la,ifirst,ilast, omega,
     &  matrix,rhs,
     &  solution,
     &  residual)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fi,la,ifirst,ilast
      double precision omega
      double precision matrix(fi:la,-1:1),rhs(fi:la)
      double precision solution(fi:la)
      double precision residual(fi:la)

      integer ic
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call compute_residual(fi,la,ifirst,ilast,matrix,rhs,solution,
     &  residual)
      do ic=ifirst,ilast
        solution(ic)=solution(ic)-omega*residual(ic)/matrix(ic,0)
      enddo
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

      subroutine gauss_seidel(fi,la,ifirst,ilast,
     &  matrix,rhs,
     &  solution)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fi,la,ifirst,ilast
      double precision matrix(fi:la,-1:1),rhs(fi:la)
      double precision solution(fi:la)

      integer ic
      double precision residual
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "in gauss_seidel, fi,la = ",fi,la
c     print *, "ifirst,ilast = ",ifirst,ilast
c     do ic=ifirst,ilast
c       print *, "matrix[",ic,"]= ",matrix(ic,-1),matrix(ic,0),
c    &    matrix(ic,1)
c     enddo
c     do ic=ifirst,ilast
c       print *, "rhs[",ic,"]= ",rhs(ic)
c     enddo
c     do ic=ifirst,ilast
c       print *, "solution[",ic,"]= ",solution(ic)
c     enddo

      ic=ifirst
      residual=matrix(ic,0)*solution(ic)
     &        +matrix(ic,1)*solution(ic+1)-rhs(ic)
      solution(ic)=solution(ic)-residual/matrix(ic,0)
      do ic=ifirst+1,ilast-1
        residual=matrix(ic,-1)*solution(ic-1)
     &          +matrix(ic, 0)*solution(ic  )
     &          +matrix(ic, 1)*solution(ic+1)
     &          -rhs(ic)
        solution(ic)=solution(ic)-residual/matrix(ic,0)
      enddo
      ic=ilast
      residual=matrix(ic,-1)*solution(ic-1)
     &        +matrix(ic, 0)*solution(ic  )-rhs(ic)
      solution(ic)=solution(ic)-residual/matrix(ic,0)

c     do ic=ifirst,ilast
c       print *, "solution[",ic,"]= ",solution(ic)
c     enddo
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

      subroutine gauss_seidel_reverse(fi,la,ifirst,ilast,
     &  matrix,rhs,
     &  solution)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fi,la,ifirst,ilast
      double precision matrix(fi:la,-1:1),rhs(fi:la)
      double precision solution(fi:la)

      integer ic
      double precision residual
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ic=ilast
      residual=matrix(ic,-1)*solution(ic-1)
     &        +matrix(ic, 0)*solution(ic  )-rhs(ic)
      solution(ic)=solution(ic)-residual/matrix(ic,0)
      do ic=ilast-1,ifirst+1,-1
        residual=matrix(ic,-1)*solution(ic-1)
     &          +matrix(ic, 0)*solution(ic  )
     &          +matrix(ic, 1)*solution(ic+1)
     &          -rhs(ic)
        solution(ic)=solution(ic)-residual/matrix(ic,0)
      enddo
      ic=ifirst
      residual=matrix(ic,0)*solution(ic)
     &        +matrix(ic,1)*solution(ic+1)-rhs(ic)
      solution(ic)=solution(ic)-residual/matrix(ic,0)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

      subroutine gauss_seidel_to_fro(fi,la,ifirst,ilast, iteration,
     &  matrix,rhs,
     &  solution)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fi,la,ifirst,ilast
      integer iteration
      double precision matrix(fi:la,-1:1),rhs(fi:la)
      double precision solution(fi:la)

      integer ic
      double precision residual
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (mod(iteration,2).eq.0) then
        ic=ifirst
        residual=matrix(ic,0)*solution(ic)
     &          +matrix(ic,1)*solution(ic+1)-rhs(ic)
        solution(ic)=solution(ic)-residual/matrix(ic,0)
        do ic=ifirst+1,ilast-1
          residual=matrix(ic,-1)*solution(ic-1)
     &            +matrix(ic, 0)*solution(ic  )
     &            +matrix(ic, 1)*solution(ic+1)
     &            -rhs(ic)
          solution(ic)=solution(ic)-residual/matrix(ic,0)
        enddo
        ic=ilast
        residual=matrix(ic,-1)*solution(ic-1)
     &          +matrix(ic, 0)*solution(ic  )-rhs(ic)
        solution(ic)=solution(ic)-residual/matrix(ic,0)
      else
        ic=ilast
        residual=matrix(ic,-1)*solution(ic-1)
     &          +matrix(ic, 0)*solution(ic  )-rhs(ic)
        solution(ic)=solution(ic)-residual/matrix(ic,0)
        do ic=ilast-1,ifirst+1,-1
          residual=matrix(ic,-1)*solution(ic-1)
     &            +matrix(ic, 0)*solution(ic  )
     &            +matrix(ic, 1)*solution(ic+1)
     &            -rhs(ic)
          solution(ic)=solution(ic)-residual/matrix(ic,0)
        enddo
        ic=ifirst
        residual=matrix(ic,0)*solution(ic)
     &          +matrix(ic,1)*solution(ic+1)-rhs(ic)
        solution(ic)=solution(ic)-residual/matrix(ic,0)
      endif
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

      subroutine gauss_seidel_red_black(fi,la,ifirst,ilast, iteration,
     &  matrix,rhs,
     &  solution)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fi,la,ifirst,ilast
      integer iteration
      double precision matrix(fi:la,-1:1),rhs(fi:la)
      double precision solution(fi:la)

      integer ic,modit
      double precision residual
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (mod(iteration,2).eq.0) then
        modit=0
        if (mod(ifirst,2).eq.modit) then
          ic=ifirst
          residual=matrix(ic,0)*solution(ic)
     &            +matrix(ic,1)*solution(ic+1)-rhs(ic)
          solution(ic)=solution(ic)-residual/matrix(ic,0)
        endif
        do ic=ifirst+1,ilast-1
          if (mod(ic,2).eq.modit) then
            residual=matrix(ic,-1)*solution(ic-1)
     &              +matrix(ic, 0)*solution(ic  )
     &              +matrix(ic, 1)*solution(ic+1)
     &              -rhs(ic)
            solution(ic)=solution(ic)-residual/matrix(ic,0)
          endif
        enddo
        if (mod(ilast,2).eq.modit) then
          ic=ilast
          residual=matrix(ic,-1)*solution(ic-1)
     &            +matrix(ic, 0)*solution(ic  )-rhs(ic)
          solution(ic)=solution(ic)-residual/matrix(ic,0)
        endif

        modit=1-modit
        if (mod(ifirst,2).eq.modit) then
          ic=ifirst
          residual=matrix(ic,0)*solution(ic)
     &            +matrix(ic,1)*solution(ic+1)-rhs(ic)
          solution(ic)=solution(ic)-residual/matrix(ic,0)
        endif
        do ic=ifirst+1,ilast-1
          if (mod(ic,2).eq.modit) then
            residual=matrix(ic,-1)*solution(ic-1)
     &              +matrix(ic, 0)*solution(ic  )
     &              +matrix(ic, 1)*solution(ic+1)
     &              -rhs(ic)
            solution(ic)=solution(ic)-residual/matrix(ic,0)
          endif
        enddo
        if (mod(ilast,2).eq.modit) then
          ic=ilast
          residual=matrix(ic,-1)*solution(ic-1)
     &            +matrix(ic, 0)*solution(ic  )-rhs(ic)
          solution(ic)=solution(ic)-residual/matrix(ic,0)
        endif
      else
        modit=0
        if (mod(ilast,2).eq.modit) then
          ic=ilast
          residual=matrix(ic,-1)*solution(ic-1)
     &            +matrix(ic, 0)*solution(ic  )-rhs(ic)
          solution(ic)=solution(ic)-residual/matrix(ic,0)
        endif
        do ic=ilast-1,ifirst+1,-1
          if (mod(ic,2).eq.modit) then
            residual=matrix(ic,-1)*solution(ic-1)
     &              +matrix(ic, 0)*solution(ic  )
     &              +matrix(ic, 1)*solution(ic+1)
     &              -rhs(ic)
            solution(ic)=solution(ic)-residual/matrix(ic,0)
          endif
        enddo
        if (mod(ifirst,2).eq.modit) then
          ic=ifirst
          residual=matrix(ic,0)*solution(ic)
     &            +matrix(ic,1)*solution(ic+1)-rhs(ic)
          solution(ic)=solution(ic)-residual/matrix(ic,0)
        endif

        modit=1-modit
        if (mod(ilast,2).eq.modit) then
          ic=ilast
          residual=matrix(ic,-1)*solution(ic-1)
     &            +matrix(ic, 0)*solution(ic  )-rhs(ic)
          solution(ic)=solution(ic)-residual/matrix(ic,0)
        endif
        do ic=ilast-1,ifirst+1,-1
          if (mod(ic,2).eq.modit) then
            residual=matrix(ic,-1)*solution(ic-1)
     &              +matrix(ic, 0)*solution(ic  )
     &              +matrix(ic, 1)*solution(ic+1)
     &              -rhs(ic)
            solution(ic)=solution(ic)-residual/matrix(ic,0)
          endif
        enddo
        if (mod(ifirst,2).eq.modit) then
          ic=ifirst
          residual=matrix(ic,0)*solution(ic)
     &            +matrix(ic,1)*solution(ic+1)-rhs(ic)
          solution(ic)=solution(ic)-residual/matrix(ic,0)
        endif
      endif
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

      subroutine sor(fi,la,ifirst,ilast, omega,
     &  matrix,rhs,
     &  solution)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fi,la,ifirst,ilast
      double precision omega
      double precision matrix(fi:la,-1:1),rhs(fi:la)
      double precision solution(fi:la)

      integer ic
      double precision residual
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ic=ifirst
      residual=matrix(ic,0)*solution(ic)
     &        +matrix(ic,1)*solution(ic+1)-rhs(ic)
      solution(ic)=solution(ic)-omega*residual/matrix(ic,0)
      do ic=ifirst+1,ilast-1
        residual=matrix(ic,-1)*solution(ic-1)
     &          +matrix(ic, 0)*solution(ic  )
     &          +matrix(ic, 1)*solution(ic+1)
     &          -rhs(ic)
        solution(ic)=solution(ic)-omega*residual/matrix(ic,0)
      enddo
      ic=ilast
      residual=matrix(ic,-1)*solution(ic-1)
     &        +matrix(ic, 0)*solution(ic  )-rhs(ic)
      solution(ic)=solution(ic)-omega*residual/matrix(ic,0)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

      subroutine algebraic_multigrid_prolongation(fic,lac,fif,laf,
     &  ifirstc,ilastc,
     &  matrixf,solutionc,
     &  solutionf)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fic,lac,fif,laf,ifirstc,ilastc
      double precision 
     &  matrixf(fif:laf,-1:1),
     &  solutionc(fic:lac)
      double precision 
     &  solutionf(fif:laf)

      integer i,ic
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "in algebraic_multigrid_prolongation, fic,lac = ",fic,lac
c     print *, "fif,laf = ",fif,laf
c     print *, "ifirstc,ilastc = ",ifirstc,ilastc
c     do i=fif,laf
c       print *, "matrixf[",i,"]= ",matrixf(i,-1),matrixf(i,0),
c    &    matrixf(i,1)
c     enddo
c     do ic=ifirstc,ilastc
c       print *, "solutionc[",ic,"]= ",solutionc(ic)
c     enddo

      do ic=ifirstc,ilastc
        i=2*ic
        solutionf(i)=solutionc(ic)
      enddo
      i=2*ifirstc-1
      solutionf(i)=-solutionf(i+1)*matrixf(i,1)/matrixf(i,0)
      do ic=ifirstc,ilastc-1
        i=2*ic+1
        solutionf(i)=-(solutionf(i-1)*matrixf(i,-1)
     &                +solutionf(i+1)*matrixf(i,1))/matrixf(i,0)
      enddo
      ic=ilastc
      i=2*ic+1
      solutionf(i)=-solutionf(i-1)*matrixf(i,-1)/matrixf(i,0)

c     do i=fif,laf
c       print *, "solutionf[",i,"]= ",solutionf(i)
c     enddo
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

      subroutine algebraic_multigrid_restriction(fic,lac,fif,laf,
     &  ifirstc,ilastc,
     &  matrixf,residualf,
     &  residualc)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fic,lac,fif,laf,ifirstc,ilastc
      double precision 
     &  matrixf(fif:laf,-1:1),
     &  residualf(fif:laf)
      double precision 
     &  residualc(fic:lac)

      integer i,ic
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      do ic=ifirstc,ilastc
        i=2*ic
        residualc(ic)=residualf(i)
     &               -residualf(i-1)*matrixf(i-1,1)/matrixf(i-1,0)
     &               -residualf(i+1)*matrixf(i+1,-1)/matrixf(i+1,0)
      enddo
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

      subroutine finite_element_prolongation(fic,lac,fif,laf,
     &  ifirstc,ilastc,
     &  solutionc,
     &  solutionf)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fic,lac,fif,laf,ifirstc,ilastc
      double precision 
     &  solutionc(fic:lac)
      double precision 
     &  solutionf(fif:laf)

      integer i,ic
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      do ic=ifirstc,ilastc
        i=2*ic
        solutionf(i)=solutionc(ic)
      enddo
      i=2*ifirstc-1
      solutionf(i)=0.5d0*solutionf(i+1)
      do ic=ifirstc,ilastc-1
        i=2*ic+1
        solutionf(i)=0.5d0*(solutionf(i-1)+solutionf(i+1))
      enddo
      ic=ilastc
      i=2*ic+1
      solutionf(i)=0.5d0*solutionf(i-1)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

      subroutine finite_element_restriction(fic,lac,fif,laf,
     &  ifirstc,ilastc,
     &  residualf,
     &  residualc)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fic,lac,fif,laf,ifirstc,ilastc
      double precision 
     &  residualf(fif:laf)
      double precision 
     &  residualc(fic:lac)

      integer i,ic
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      do ic=ifirstc,ilastc
        i=2*ic
        residualc(ic)=
     &    residualf(i)+0.5d0*(residualf(i-1)+residualf(i+1))
      enddo
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 
