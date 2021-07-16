c "$Header: /home/faculty/johnt/cvs/memdebug/m4amr.sed,v 1.11 2003/06/15 13:50:01 johnt Exp $"






c assume Dirichlet boundary conditions in all methods ==>
c   no change to solution at boundary

      subroutine setup_heat(fi_0,fi_1,la_0,la_1,
     &  ifirst_0,ifirst_1,ilast_0,ilast_1,decay_0,
     &  diffusion,true_solution,
     &  matrix,solution)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c***********************************************************************
c  Copyright 2006 John A. Trangenstein
c
c  This software is made available for research and instructional use 
c  only. 
c  You may copy and use this software without charge for these 
c  non-commercial purposes, provided that the copyright notice and 
c  associated text is reproduced on all copies.  
c  For all other uses (including distribution of modified versions), 
c  please contact the author at
c    John A. Trangenstein
c    Department of Mathematics
c    Duke University
c    Durham, NC 27708-0320
c    USA
c  or
c    johnt@math.duke.edu
c  
c  This software is made available "as is" without any assurance that it
c  is completely correct, or that it will work for your purposes.  
c  Use the software at your own risk.
c***********************************************************************
      double precision roundoff,small,huge,undefind
      common/machine/roundoff,small,huge,undefind

      double precision e,pi,rt2
      double precision half,one,two,zero
      integer debug_on
      common/const/e,pi,rt2, half,one,two,zero,debug_on
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fi_0,fi_1,la_0,la_1,ifirst_0,ifirst_1,ilast_0,ilast_1
      double precision decay_0,decay_1
      double precision 
     &  diffusion(fi_0:la_0-1,fi_1:la_1-1),
     &  true_solution(fi_0:la_0,fi_1:la_1)
      double precision 
     &  matrix(fi_0:la_0,fi_1:la_1,-1:1,-1:1),
     &  solution(fi_0:la_0,fi_1:la_1)

      integer ic_0,ic_1
      double precision dijm,dijp,dimj,dipj,dt,dx_0,dx_1
      integer i,j
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "in setup_heat"

      dx_0=1.d0/dble(la_0-fi_0)
      dx_1=1.d0/dble(la_1-fi_1)
      dt=decay_0*dx_0**2
      decay_1=decay_0*(dx_0/dx_1)**2

      do j=-1,1
        do i=-1,1
          do ic_1=ifirst_1,ilast_1
            do ic_0=ifirst_0,ilast_0
              matrix(ic_0,ic_1,i,j)=0.d0
            enddo
          enddo
        enddo
      enddo

c     finite difference laplacian; assumes Dirichlet bc's all around
      do ic_1=ifirst_1,ilast_1
        do ic_0=ifirst_0,ilast_0
          dipj=0.5d0*(diffusion(ic_0  ,ic_1-1)+diffusion(ic_0  ,ic_1  ))
          dimj=0.5d0*(diffusion(ic_0-1,ic_1-1)+diffusion(ic_0-1,ic_1  ))
          dijp=0.5d0*(diffusion(ic_0-1,ic_1  )+diffusion(ic_0  ,ic_1  ))
          dijm=0.5d0*(diffusion(ic_0-1,ic_1-1)+diffusion(ic_0  ,ic_1-1))
          matrix(ic_0,ic_1,1,0)=-dt*dipj*dx_1/dx_0
          matrix(ic_0,ic_1,-1,0)=-dt*dimj*dx_1/dx_0
          matrix(ic_0,ic_1,0,1)=-dt*dijp*dx_0/dx_1
          matrix(ic_0,ic_1,0,-1)=-dt*dijm*dx_0/dx_1
          matrix(ic_0,ic_1,0,0)=dx_0*dx_1
     &      -matrix(ic_0,ic_1,1,0)-matrix(ic_0,ic_1,-1,0)
     &      -matrix(ic_0,ic_1,0,1)-matrix(ic_0,ic_1,0,-1)
        enddo
      enddo

c     set Dirichlet boundary values for graphics
      do ic_1=ifirst_1,ilast_1
        solution(fi_0,ic_1)=true_solution(fi_0,ic_1)
        solution(la_0,ic_1)=true_solution(la_0,ic_1)
      enddo
      do ic_0=fi_0,la_0
        solution(ic_0,fi_1)=true_solution(ic_0,fi_1)
        solution(ic_0,la_1)=true_solution(ic_0,la_1)
      enddo

c     print *, "leaving setup_heat"
c     print *, "solution = "
c     call rbugcell(fi_0,fi_1,la_0,la_1,
c    &  ifirst_0,ifirst_1,ilast_0,ilast_1,solution)
c     do j=-1,1
c       do i=-1,1
c         print *, "matrix[",i,",",j,"] = "
c         call rbugcell(fi_0,fi_1,la_0,la_1,
c    &      ifirst_0,ifirst_1,ilast_0,ilast_1,matrix(fi_0,fi_1,i,j))
c       enddo
c     enddo
      call flush(6)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

      subroutine setup_laplace(fi_0,fi_1,la_0,la_1,
     &  ifirst_0,ifirst_1,ilast_0,ilast_1,
     &  diffusion,true_solution,
     &  matrix,solution)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c***********************************************************************
c  Copyright 2006 John A. Trangenstein
c
c  This software is made available for research and instructional use 
c  only. 
c  You may copy and use this software without charge for these 
c  non-commercial purposes, provided that the copyright notice and 
c  associated text is reproduced on all copies.  
c  For all other uses (including distribution of modified versions), 
c  please contact the author at
c    John A. Trangenstein
c    Department of Mathematics
c    Duke University
c    Durham, NC 27708-0320
c    USA
c  or
c    johnt@math.duke.edu
c  
c  This software is made available "as is" without any assurance that it
c  is completely correct, or that it will work for your purposes.  
c  Use the software at your own risk.
c***********************************************************************
      double precision roundoff,small,huge,undefind
      common/machine/roundoff,small,huge,undefind

      double precision e,pi,rt2
      double precision half,one,two,zero
      integer debug_on
      common/const/e,pi,rt2, half,one,two,zero,debug_on
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fi_0,fi_1,la_0,la_1,ifirst_0,ifirst_1,ilast_0,ilast_1
      double precision 
     &  diffusion(fi_0:la_0-1,fi_1:la_1-1),
     &  true_solution(fi_0:la_0,fi_1:la_1)
      double precision 
     &  matrix(fi_0:la_0,fi_1:la_1,-1:1,-1:1),
     &  solution(fi_0:la_0,fi_1:la_1)

      integer ic_0,ic_1
      double precision dijm,dijp,dimj,dipj,dx_0,dx_1
c     double precision a,ad,h,hd,s,sd
      integer i,j
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "in setup_laplace"
c     print *, "diffusion = "
c     call rbugcell(fi_0,fi_1,la_0-1,la_1-1,
c    &  fi_0,fi_1,la_0-1,la_1-1,diffusion)
c     call flush(6)

      do j=-1,1
        do i=-1,1
          do ic_1=ifirst_1,ilast_1
            do ic_0=ifirst_0,ilast_0
              matrix(ic_0,ic_1,i,j)=0.d0
            enddo
          enddo
        enddo
      enddo

      dx_0=1.d0/dble(la_0-fi_0)
      dx_1=1.d0/dble(la_1-fi_1)

c     finite difference laplacian; assumes Dirichlet bc's all around
      do ic_1=ifirst_1,ilast_1
        do ic_0=ifirst_0,ilast_0
          dipj=0.5d0*(diffusion(ic_0  ,ic_1-1)+diffusion(ic_0  ,ic_1  ))
          dimj=0.5d0*(diffusion(ic_0-1,ic_1-1)+diffusion(ic_0-1,ic_1  ))
          dijp=0.5d0*(diffusion(ic_0-1,ic_1  )+diffusion(ic_0  ,ic_1  ))
          dijm=0.5d0*(diffusion(ic_0-1,ic_1-1)+diffusion(ic_0  ,ic_1-1))
          matrix(ic_0,ic_1,1,0)=-dipj*dx_1/dx_0
          matrix(ic_0,ic_1,-1,0)=-dimj*dx_1/dx_0
          matrix(ic_0,ic_1,0,1)=-dijp*dx_0/dx_1
          matrix(ic_0,ic_1,0,-1)=-dijm*dx_0/dx_1
          matrix(ic_0,ic_1,0,0)=
     &      -matrix(ic_0,ic_1,1,0)-matrix(ic_0,ic_1,-1,0)
     &      -matrix(ic_0,ic_1,0,1)-matrix(ic_0,ic_1,0,-1)
        enddo
      enddo

c     set Dirichlet boundary values for graphics
      do ic_1=ifirst_1,ilast_1
        solution(fi_0,ic_1)=true_solution(fi_0,ic_1)
        solution(la_0,ic_1)=true_solution(la_0,ic_1)
      enddo
      do ic_0=fi_0,la_0
        solution(ic_0,fi_1)=true_solution(ic_0,fi_1)
        solution(ic_0,la_1)=true_solution(ic_0,la_1)
      enddo

c     print *, "leaving setup_laplace"
c     print *, "solution = "
c     call rbugcell(fi_0,fi_1,la_0,la_1,
c    &  ifirst_0,ifirst_1,ilast_0,ilast_1,solution)
c     do j=-1,1
c       do i=-1,1
c         print *, "matrix[",i,",",j,"] = "
c         call rbugcell(fi_0,fi_1,la_0,la_1,
c    &      ifirst_0,ifirst_1,ilast_0,ilast_1,matrix(fi_0,fi_1,i,j))
c       enddo
c     enddo
c     call flush(6)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

      subroutine setup_convection_diffusion(fi_0,fi_1,la_0,la_1,
     &  ifirst_0,ifirst_1,ilast_0,ilast_1,decay,peclet,
     &  true_solution,
     &  matrix,solution)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c***********************************************************************
c  Copyright 2006 John A. Trangenstein
c
c  This software is made available for research and instructional use 
c  only. 
c  You may copy and use this software without charge for these 
c  non-commercial purposes, provided that the copyright notice and 
c  associated text is reproduced on all copies.  
c  For all other uses (including distribution of modified versions), 
c  please contact the author at
c    John A. Trangenstein
c    Department of Mathematics
c    Duke University
c    Durham, NC 27708-0320
c    USA
c  or
c    johnt@math.duke.edu
c  
c  This software is made available "as is" without any assurance that it
c  is completely correct, or that it will work for your purposes.  
c  Use the software at your own risk.
c***********************************************************************
      double precision roundoff,small,huge,undefind
      common/machine/roundoff,small,huge,undefind

      double precision e,pi,rt2
      double precision half,one,two,zero
      integer debug_on
      common/const/e,pi,rt2, half,one,two,zero,debug_on
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fi_0,fi_1,la_0,la_1,ifirst_0,ifirst_1,ilast_0,ilast_1
      double precision decay,peclet
      double precision 
     &  true_solution(fi_0:la_0,fi_1:la_1)
      double precision 
     &  matrix(fi_0:la_0,fi_1:la_1,-1:1,-1:1),
     &  solution(fi_0:la_0,fi_1:la_1)

      integer ic_0,ic_1
      double precision diag,dt,subdiagi,subdiagj,supdiagi,supdiagj,
     &  dx_0,dx_1,velocity
c     integer i,j
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "in setup_convection_diffusion"
c     print *, "decay,peclet = ",decay,peclet
c     print *, "true_solution = "
c     call rbugcell(fi_0,fi_1,la_0,la_1,
c    &  ifirst_0,ifirst_1,ilast_0,ilast_1,true_solution)

      dx_0=1.d0/dble(la_0-fi_0)
      dx_1=1.d0/dble(la_1-fi_1)
      dt=decay*dx_0**2
      velocity=peclet*dx_0

      supdiagi=-dt*dx_1/dx_0
      supdiagj=-dt*dx_0/dx_1
      subdiagi=-dt*(dx_1/dx_0+velocity*dx_1)
      subdiagj=-dt*(dx_0/dx_1+velocity*dx_0)
      diag=dx_0*dx_1+dt*(2.*(dx_1/dx_0+dx_0/dx_1)+velocity*(dx_1+dx_0))
      do ic_1=ifirst_1,ilast_1
        do ic_0=ifirst_0,ilast_0
c         backward euler in time:
          matrix(ic_0,ic_1,1,0)=supdiagi
          matrix(ic_0,ic_1,-1,0)=subdiagi
          matrix(ic_0,ic_1,0,1)=supdiagj
          matrix(ic_0,ic_1,0,-1)=subdiagj
          matrix(ic_0,ic_1,0,0)=diag
          matrix(ic_0,ic_1,-1,1)=0.d0
          matrix(ic_0,ic_1, 1,1)=0.d0
          matrix(ic_0,ic_1,-1,-1)=0.d0
          matrix(ic_0,ic_1, 1,-1)=0.d0
        enddo
c       dirichlet boundary condition:
        solution(fi_0,ic_1)=true_solution(fi_0,ic_1)
        solution(la_0,ic_1)=true_solution(la_0,ic_1)
      enddo
c     dirichlet boundary condition:
      do ic_0=fi_0,la_0
        solution(ic_0,fi_1)=true_solution(ic_0,fi_1)
        solution(ic_0,la_1)=true_solution(ic_0,la_1)
      enddo

c     print *, "leaving setup_convection_diffusion"
c     print *, "solution = "
c     call rbugcell(fi_0,fi_1,la_0,la_1,
c    &  ifirst_0,ifirst_1,ilast_0,ilast_1,solution)
c     do j=-1,1
c       do i=-1,1
c         print *, "matrix[",i,",",j,"] = "
c         call rbugcell(fi_0,fi_1,la_0,la_1,
c    &      ifirst_0,ifirst_1,ilast_0,ilast_1,matrix(fi_0,fi_1,i,j))
c       enddo
c     enddo
      call flush(6)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

      subroutine setup_dlap_matrix(fi_0,fi_1,la_0,la_1,
     &  ifirst_0,ifirst_1,ilast_0,ilast_1,nelt,
     &  matrix,
     &  dlap_a,dlap_ia,dlap_ja)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c***********************************************************************
c  Copyright 2006 John A. Trangenstein
c
c  This software is made available for research and instructional use 
c  only. 
c  You may copy and use this software without charge for these 
c  non-commercial purposes, provided that the copyright notice and 
c  associated text is reproduced on all copies.  
c  For all other uses (including distribution of modified versions), 
c  please contact the author at
c    John A. Trangenstein
c    Department of Mathematics
c    Duke University
c    Durham, NC 27708-0320
c    USA
c  or
c    johnt@math.duke.edu
c  
c  This software is made available "as is" without any assurance that it
c  is completely correct, or that it will work for your purposes.  
c  Use the software at your own risk.
c***********************************************************************
      double precision roundoff,small,huge,undefind
      common/machine/roundoff,small,huge,undefind

      double precision e,pi,rt2
      double precision half,one,two,zero
      integer debug_on
      common/const/e,pi,rt2, half,one,two,zero,debug_on
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer 
     &  fi_0,fi_1,la_0,la_1,
     &  ifirst_0,ifirst_1,ilast_0,ilast_1,
     &  nelt
      double precision 
     &  matrix(fi_0:la_0,fi_1:la_1,-1:1,-1:1)
      double precision 
     &  dlap_a(1)
      integer
     &  dlap_ia(1),
     &  dlap_ja(1)

      integer ic_0,ic_1,n_0,n_1,row
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "in setup_dlap_matrix"
c     print *, "fi = ",fi_0,fi_1
c     print *, "la = ",la_0,la_1
c     print *, "ifirst = ",ifirst_0,ifirst_1
c     print *, "ilast = ",ilast_0,ilast_1

      nelt=1
      n_0=ilast_0-ifirst_0+1
      n_1=ilast_1-ifirst_1+1
c     print *, "n = ",n_0,n_1

      row=1
      do ic_1=ifirst_1,ilast_1
        do ic_0=ifirst_0,ilast_0
c         print *
c         print *, "ic = ",ic_0,ic_1
          if (ic_1.gt.ifirst_1) then
            if (ic_0.gt.ifirst_0) then
              dlap_a(nelt)=matrix(ic_0,ic_1,-1,-1)
              dlap_ia(nelt)=row
              dlap_ja(nelt)=row-1-n_0
c             print *, "using -1,-1: a[",dlap_ia(nelt),",",
c    &          dlap_ja(nelt),"] = ",dlap_a(nelt)
              nelt=nelt+1
            endif
            dlap_a(nelt)=matrix(ic_0,ic_1,0,-1)
            dlap_ia(nelt)=row
            dlap_ja(nelt)=row-n_0
c           print *, "using 0,-1: a[",dlap_ia(nelt),",",dlap_ja(nelt),
c    &          "] = ",dlap_a(nelt)
            nelt=nelt+1
            if (ic_0.lt.ilast_0) then
              dlap_a(nelt)=matrix(ic_0,ic_1,1,-1)
              dlap_ia(nelt)=row
              dlap_ja(nelt)=row+1-n_0
c             print *, "using 1,-1: a[",dlap_ia(nelt),",",
c    &          dlap_ja(nelt),"] = ",dlap_a(nelt)
              nelt=nelt+1
            endif
          endif
          if (ic_0.gt.ifirst_0) then
            dlap_a(nelt)=matrix(ic_0,ic_1,-1,0)
            dlap_ia(nelt)=row
            dlap_ja(nelt)=row-1
c           print *, "using -1,0: a[",dlap_ia(nelt),",",dlap_ja(nelt),
c    &          "] = ",dlap_a(nelt)
            nelt=nelt+1
          endif
          dlap_a(nelt)=matrix(ic_0,ic_1,0,0)
          dlap_ia(nelt)=row
          dlap_ja(nelt)=row
c         print *, "using 0,0: a[",dlap_ia(nelt),",",dlap_ja(nelt),
c    &          "] = ",dlap_a(nelt)," = matrix[",ic_0,",",ic_1,",0,0]"
          nelt=nelt+1
          if (ic_0.lt.ilast_0) then
            dlap_a(nelt)=matrix(ic_0,ic_1,1,0)
            dlap_ia(nelt)=row
            dlap_ja(nelt)=row+1
c           print *, "using 1,0: a[",dlap_ia(nelt),",",dlap_ja(nelt),
c    &          "] = ",dlap_a(nelt)
            nelt=nelt+1
          endif
          if (ic_1.lt.ilast_1) then
            if (ic_0.gt.ifirst_0) then
              dlap_a(nelt)=matrix(ic_0,ic_1,-1,1)
              dlap_ia(nelt)=row
              dlap_ja(nelt)=row-1+n_0
c             print *, "using -1,1: a[",dlap_ia(nelt),",",
c    &          dlap_ja(nelt),"] = ",dlap_a(nelt)
              nelt=nelt+1
            endif
            dlap_a(nelt)=matrix(ic_0,ic_1,0,1)
            dlap_ia(nelt)=row
            dlap_ja(nelt)=row+n_0
c           print *, "using 0,1: a[",dlap_ia(nelt),",",dlap_ja(nelt),
c    &          "] = ",dlap_a(nelt)
            nelt=nelt+1
            if (ic_0.lt.ilast_0) then
              dlap_a(nelt)=matrix(ic_0,ic_1,1,1)
              dlap_ia(nelt)=row
              dlap_ja(nelt)=row+1+n_0
c             print *, "using 1,1: a[",dlap_ia(nelt),",",
c    &          dlap_ja(nelt),"] = ",dlap_a(nelt)
              nelt=nelt+1
            endif
          endif
          row=row+1
        enddo
      enddo
      nelt=nelt-1

c     print *, "leaving setup_dlap_matrix, nelt = ",nelt
c     do i=1,nelt
c       print *, "i,a,ia,ja = ",i,dlap_a(i),dlap_ia(i),dlap_ja(i)
c     enddo
c     call flush(6)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

      subroutine setup_dlap_system(fi_0,fi_1,la_0,la_1,
     &  ifirst_0,ifirst_1,ilast_0,ilast_1,
     &  x,
     &  dlap_x)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c***********************************************************************
c  Copyright 2006 John A. Trangenstein
c
c  This software is made available for research and instructional use 
c  only. 
c  You may copy and use this software without charge for these 
c  non-commercial purposes, provided that the copyright notice and 
c  associated text is reproduced on all copies.  
c  For all other uses (including distribution of modified versions), 
c  please contact the author at
c    John A. Trangenstein
c    Department of Mathematics
c    Duke University
c    Durham, NC 27708-0320
c    USA
c  or
c    johnt@math.duke.edu
c  
c  This software is made available "as is" without any assurance that it
c  is completely correct, or that it will work for your purposes.  
c  Use the software at your own risk.
c***********************************************************************
      double precision roundoff,small,huge,undefind
      common/machine/roundoff,small,huge,undefind

      double precision e,pi,rt2
      double precision half,one,two,zero
      integer debug_on
      common/const/e,pi,rt2, half,one,two,zero,debug_on
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer 
     &  fi_0,fi_1,la_0,la_1,
     &  ifirst_0,ifirst_1,ilast_0,ilast_1
      double precision 
     &  x(fi_0:la_0,fi_1:la_1)
      double precision 
     &  dlap_x(1)

      integer ic_0,ic_1,k
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "in setup_dlap_system"
c     print *, "x = "
c     call rbugcell(fi_0,fi_1,la_0,la_1,
c    &  ifirst_0,ifirst_1,ilast_0,ilast_1,x)

      k=1
      do ic_1=ifirst_1,ilast_1
        do ic_0=ifirst_0,ilast_0
          dlap_x(k)=x(ic_0,ic_1)
          k=k+1
        enddo
      enddo

c     print *, "leaving setup_dlap_system"
c     do i=1,k-1
c       print *, "dlap_x[",i,"] = ",dlap_x(i)
c     enddo
c     call flush(6)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

      subroutine copy_dlap_solution(fi_0,fi_1,la_0,la_1,
     &  ifirst_0,ifirst_1,ilast_0,ilast_1,
     &  dlap_x,
     &  x)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c***********************************************************************
c  Copyright 2006 John A. Trangenstein
c
c  This software is made available for research and instructional use 
c  only. 
c  You may copy and use this software without charge for these 
c  non-commercial purposes, provided that the copyright notice and 
c  associated text is reproduced on all copies.  
c  For all other uses (including distribution of modified versions), 
c  please contact the author at
c    John A. Trangenstein
c    Department of Mathematics
c    Duke University
c    Durham, NC 27708-0320
c    USA
c  or
c    johnt@math.duke.edu
c  
c  This software is made available "as is" without any assurance that it
c  is completely correct, or that it will work for your purposes.  
c  Use the software at your own risk.
c***********************************************************************
      double precision roundoff,small,huge,undefind
      common/machine/roundoff,small,huge,undefind

      double precision e,pi,rt2
      double precision half,one,two,zero
      integer debug_on
      common/const/e,pi,rt2, half,one,two,zero,debug_on
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer 
     &  fi_0,fi_1,la_0,la_1,
     &  ifirst_0,ifirst_1,ilast_0,ilast_1
      double precision 
     &  dlap_x(1)
      double precision 
     &  x(fi_0:la_0,fi_1:la_1)

      integer ic_0,ic_1,k
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "in copy_dlap_solution"

      k=1
      do ic_1=ifirst_1,ilast_1
        do ic_0=ifirst_0,ilast_0
          x(ic_0,ic_1)=dlap_x(k)
          k=k+1
        enddo
      enddo

c     print *, "leaving copy_dlap_solution"
c     do i=1,k-1
c       print *, "dlap_x[",i,"] = ",dlap_x(i)
c     enddo
c     print *, "x = "
c     call rbugcell(fi_0,fi_1,la_0,la_1,
c    &  ifirst_0,ifirst_1,ilast_0,ilast_1,x)
c     call flush(6)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

      function max_residual(fi_0,fi_1,la_0,la_1,
     &  ifirst_0,ifirst_1,ilast_0,ilast_1,
     &  residual)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c***********************************************************************
c  Copyright 2006 John A. Trangenstein
c
c  This software is made available for research and instructional use 
c  only. 
c  You may copy and use this software without charge for these 
c  non-commercial purposes, provided that the copyright notice and 
c  associated text is reproduced on all copies.  
c  For all other uses (including distribution of modified versions), 
c  please contact the author at
c    John A. Trangenstein
c    Department of Mathematics
c    Duke University
c    Durham, NC 27708-0320
c    USA
c  or
c    johnt@math.duke.edu
c  
c  This software is made available "as is" without any assurance that it
c  is completely correct, or that it will work for your purposes.  
c  Use the software at your own risk.
c***********************************************************************
      double precision roundoff,small,huge,undefind
      common/machine/roundoff,small,huge,undefind

      double precision e,pi,rt2
      double precision half,one,two,zero
      integer debug_on
      common/const/e,pi,rt2, half,one,two,zero,debug_on
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fi_0,fi_1,la_0,la_1,ifirst_0,ifirst_1,ilast_0,ilast_1
      double precision max_residual, 
     &  residual(fi_0:la_0,fi_1:la_1)

      integer ic_0,ic_1
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      max_residual=0.d0
      do ic_1=ifirst_1,ilast_1
        do ic_0=ifirst_0,ilast_0
          max_residual=max(max_residual,abs(residual(ic_0,ic_1)))
        enddo
      enddo
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

      subroutine compute_residual(fi_0,fi_1,la_0,la_1,
     &  ifirst_0,ifirst_1,ilast_0,ilast_1, residual_max,
     &  matrix,rhs,solution,
     &  residual)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fi_0,fi_1,la_0,la_1,ifirst_0,ifirst_1,ilast_0,ilast_1
      double precision residual_max 
      double precision 
     &  matrix(fi_0:la_0,fi_1:la_1,-1:1,-1:1),
     &  rhs(fi_0:la_0,fi_1:la_1),
     &  solution(fi_0:la_0,fi_1:la_1)
      double precision residual(fi_0:la_0,fi_1:la_1)

      integer ic_0,ic_1
c     integer i,j
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "in compute_residual"
c     print *, "fi = ",fi_0,fi_1
c     print *, "la = ",la_0,la_1
c     print *, "ifirst = ",ifirst_0,ifirst_1
c     print *, "ilast = ",ilast_0,ilast_1
c     print *, "rhs = "
c     call rbugcell(fi_0,fi_1,la_0,la_1,
c    &  ifirst_0,ifirst_1,ilast_0,ilast_1,rhs)
c     print *, "solution = "
c     call rbugcell(fi_0,fi_1,la_0,la_1,
c    &  ifirst_0,ifirst_1,ilast_0,ilast_1,solution)
c     do j=-1,1
c       do i=-1,1
c         print *, "matrix[",i,",",j,"] = "
c         call rbugcell(fi_0,fi_1,la_0,la_1,
c    &      ifirst_0,ifirst_1,ilast_0,ilast_1,matrix(fi_0,fi_1,i,j))
c       enddo
c     enddo
c     call flush(6)

      residual_max=0.d0
      do ic_1=ifirst_1,ilast_1
        do ic_0=ifirst_0,ilast_0
          residual(ic_0,ic_1)=
     &       matrix(ic_0,ic_1,-1,-1)*solution(ic_0-1,ic_1-1)
     &      +matrix(ic_0,ic_1, 0,-1)*solution(ic_0  ,ic_1-1)
     &      +matrix(ic_0,ic_1, 1,-1)*solution(ic_0+1,ic_1-1)
     &      +matrix(ic_0,ic_1,-1, 0)*solution(ic_0-1,ic_1  )
     &      +matrix(ic_0,ic_1, 0, 0)*solution(ic_0  ,ic_1  )
     &      +matrix(ic_0,ic_1, 1, 0)*solution(ic_0+1,ic_1  )
     &      +matrix(ic_0,ic_1,-1, 1)*solution(ic_0-1,ic_1+1)
     &      +matrix(ic_0,ic_1, 0, 1)*solution(ic_0  ,ic_1+1)
     &      +matrix(ic_0,ic_1, 1, 1)*solution(ic_0+1,ic_1+1)
     &      -rhs(ic_0,ic_1)
          residual_max=max(residual_max,abs(residual(ic_0,ic_1)))
        enddo
      enddo

c     print *, "leaving compute_residual"
c     print *, "residual = "
c     call rbugcell(fi_0,fi_1,la_0,la_1,
c    &  ifirst_0,ifirst_1,ilast_0,ilast_1,residual)
c     call flush(6)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

      subroutine richardson(fi_0,fi_1,la_0,la_1,
     &  ifirst_0,ifirst_1,ilast_0,ilast_1, mu,residual_max,
     &  matrix,rhs,
     &  solution,
     &  residual)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fi_0,fi_1,la_0,la_1,ifirst_0,ifirst_1,ilast_0,ilast_1
      double precision mu,residual_max
      double precision 
     &  matrix(fi_0:la_0,fi_1:la_1,-1:1,-1:1),
     &  rhs(fi_0:la_0,fi_1:la_1)
      double precision solution(fi_0:la_0,fi_1:la_1)
      double precision residual(fi_0:la_0,fi_1:la_1)

      integer ic_0,ic_1
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "in richardson"
c     print *, "rhs = "
c     call rbugcell(fi_0,fi_1,la_0,la_1,
c    &  ifirst_0,ifirst_1,ilast_0,ilast_1,rhs)
c     print *, "solution = "
c     call rbugcell(fi_0,fi_1,la_0,la_1,
c    &  ifirst_0,ifirst_1,ilast_0,ilast_1,solution)
c     do j=-1,1
c       do i=-1,1
c         print *, "matrix[",i,",",j,"] = "
c         call rbugcell(fi_0,fi_1,la_0,la_1,
c    &      ifirst_0,ifirst_1,ilast_0,ilast_1,matrix(fi_0,fi_1,i,j))
c       enddo
c     enddo
c     call flush(6)

      call compute_residual(fi_0,fi_1,la_0,la_1,
     &  ifirst_0,ifirst_1,ilast_0,ilast_1, residual_max,
     &  matrix,rhs,solution,
     &  residual)

      do ic_1=ifirst_1,ilast_1
        do ic_0=ifirst_0,ilast_0
          solution(ic_0,ic_1)=solution(ic_0,ic_1)-residual(ic_0,ic_1)/mu
        enddo
      enddo

c     print *, "leaving richardson"
c     print *, "solution = "
c     call rbugcell(fi_0,fi_1,la_0,la_1,
c    &  ifirst_0,ifirst_1,ilast_0,ilast_1,solution)
c     call flush(6)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

      subroutine jacobi(fi_0,fi_1,la_0,la_1,
     &  ifirst_0,ifirst_1,ilast_0,ilast_1,
     &  matrix,rhs,
     &  solution,
     &  residual)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fi_0,fi_1,la_0,la_1,ifirst_0,ifirst_1,ilast_0,ilast_1
      double precision 
     &  matrix(fi_0:la_0,fi_1:la_1,-1:1,-1:1),
     &  rhs(fi_0:la_0,fi_1:la_1)
      double precision solution(fi_0:la_0,fi_1:la_1)
      double precision residual(fi_0:la_0,fi_1:la_1)

      integer ic_0,ic_1
      double precision residual_max
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call compute_residual(fi_0,fi_1,la_0,la_1,
     &  ifirst_0,ifirst_1,ilast_0,ilast_1, residual_max,
     &  matrix,rhs,solution,
     &  residual)

      do ic_1=ifirst_1,ilast_1
        do ic_0=ifirst_0,ilast_0
          solution(ic_0,ic_1)=solution(ic_0,ic_1)
     &                       -residual(ic_0,ic_1)/matrix(ic_0,ic_1,0,0)
        enddo
      enddo
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

      subroutine jacobi_omega(fi_0,fi_1,la_0,la_1,
     &  ifirst_0,ifirst_1,ilast_0,ilast_1, omega,residual_max,
     &  matrix,rhs,
     &  solution,
     &  residual)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fi_0,fi_1,la_0,la_1,ifirst_0,ifirst_1,ilast_0,ilast_1
      double precision omega,residual_max
      double precision 
     &  matrix(fi_0:la_0,fi_1:la_1,-1:1,-1:1),
     &  rhs(fi_0:la_0,fi_1:la_1)
      double precision solution(fi_0:la_0,fi_1:la_1)
      double precision residual(fi_0:la_0,fi_1:la_1)

      integer ic_0,ic_1
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call compute_residual(fi_0,fi_1,la_0,la_1,
     &  ifirst_0,ifirst_1,ilast_0,ilast_1, residual_max,
     &  matrix,rhs,solution,
     &  residual)

      do ic_1=ifirst_1,ilast_1
        do ic_0=ifirst_0,ilast_0
          solution(ic_0,ic_1)=solution(ic_0,ic_1)
     &                  -omega*residual(ic_0,ic_1)/matrix(ic_0,ic_1,0,0)
        enddo
      enddo
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

      subroutine gauss_seidel_to_fro(fi_0,fi_1,la_0,la_1,
     &  ifirst_0,ifirst_1,ilast_0,ilast_1, iteration,residual_max,
     &  matrix,rhs,
     &  solution)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fi_0,fi_1,la_0,la_1,ifirst_0,ifirst_1,ilast_0,ilast_1
      integer iteration
      double precision residual_max 
      double precision 
     &  matrix(fi_0:la_0,fi_1:la_1,-1:1,-1:1),
     &  rhs(fi_0:la_0,fi_1:la_1)
      double precision solution(fi_0:la_0,fi_1:la_1)

      integer ic_0,ic_1
      double precision residual
c     integer i,j
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "in gauss_seidel_to_fro, iteration = ",iteration
c     print *, "rhs"
c     call rbugcell(fi_0,fi_1,la_0,la_1,
c    &  ifirst_0,ifirst_1,ilast_0,ilast_1,rhs)
c     print *, "solution"
c     call rbugcell(fi_0,fi_1,la_0,la_1,
c    &  ifirst_0,ifirst_1,ilast_0,ilast_1,solution)
c     do j=-1,1
c       do i=-1,1
c         print *, "matrix(",i,",",j,")"
c         call rbugcell(fi_0,fi_1,la_0,la_1,
c    &      ifirst_0,ifirst_1,ilast_0,ilast_1,matrix(fi_0,fi_1,i,j))
c       enddo
c     enddo
      residual_max=0.d0
      if (mod(iteration,2).eq.0) then
        do ic_1=ifirst_1,ilast_1
          do ic_0=ifirst_0,ilast_0
            residual=
     &       matrix(ic_0,ic_1,-1,-1)*solution(ic_0-1,ic_1-1)
     &      +matrix(ic_0,ic_1, 0,-1)*solution(ic_0  ,ic_1-1)
     &      +matrix(ic_0,ic_1, 1,-1)*solution(ic_0+1,ic_1-1)
     &      +matrix(ic_0,ic_1,-1, 0)*solution(ic_0-1,ic_1  )
     &      +matrix(ic_0,ic_1, 0, 0)*solution(ic_0  ,ic_1  )
     &      +matrix(ic_0,ic_1, 1, 0)*solution(ic_0+1,ic_1  )
     &      +matrix(ic_0,ic_1,-1, 1)*solution(ic_0-1,ic_1+1)
     &      +matrix(ic_0,ic_1, 0, 1)*solution(ic_0  ,ic_1+1)
     &      +matrix(ic_0,ic_1, 1, 1)*solution(ic_0+1,ic_1+1)
     &      -rhs(ic_0,ic_1)
            residual_max=max(residual_max,abs(residual))
            solution(ic_0,ic_1)=solution(ic_0,ic_1)
     &                         -residual/matrix(ic_0,ic_1,0,0)
          enddo
        enddo
      else
        do ic_1=ilast_1,ifirst_1,-1
          do ic_0=ilast_0,ifirst_0,-1
            residual=
     &       matrix(ic_0,ic_1,-1,-1)*solution(ic_0-1,ic_1-1)
     &      +matrix(ic_0,ic_1, 0,-1)*solution(ic_0  ,ic_1-1)
     &      +matrix(ic_0,ic_1, 1,-1)*solution(ic_0+1,ic_1-1)
     &      +matrix(ic_0,ic_1,-1, 0)*solution(ic_0-1,ic_1  )
     &      +matrix(ic_0,ic_1, 0, 0)*solution(ic_0  ,ic_1  )
     &      +matrix(ic_0,ic_1, 1, 0)*solution(ic_0+1,ic_1  )
     &      +matrix(ic_0,ic_1,-1, 1)*solution(ic_0-1,ic_1+1)
     &      +matrix(ic_0,ic_1, 0, 1)*solution(ic_0  ,ic_1+1)
     &      +matrix(ic_0,ic_1, 1, 1)*solution(ic_0+1,ic_1+1)
     &      -rhs(ic_0,ic_1)
            residual_max=max(residual_max,abs(residual))
            solution(ic_0,ic_1)=solution(ic_0,ic_1)
     &                         -residual/matrix(ic_0,ic_1,0,0)
          enddo
        enddo
      endif
c     print *, "leaving gauss_seidel_to_fro"
c     print *, "solution"
c     call rbugcell(fi_0,fi_1,la_0,la_1,
c    &  ifirst_0,ifirst_1,ilast_0,ilast_1,solution)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

      subroutine gauss_seidel_red_black(fi_0,fi_1,la_0,la_1,
     &  ifirst_0,ifirst_1,ilast_0,ilast_1, iteration,residual_max,
     &  matrix,rhs,
     &  solution)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fi_0,fi_1,la_0,la_1,ifirst_0,ifirst_1,ilast_0,ilast_1
      integer iteration
      double precision residual_max 
      double precision 
     &  matrix(fi_0:la_0,fi_1:la_1,-1:1,-1:1),
     &  rhs(fi_0:la_0,fi_1:la_1)
      double precision solution(fi_0:la_0,fi_1:la_1)

      integer ic_0,ic_1,modit,istart_0
      double precision residual
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "in gauss_seidel_red_black"
      residual_max=0.d0
      if (mod(iteration,2).eq.0) then
c       print *, "even iteration"
        modit=0
        do ic_1=ifirst_1,ilast_1
          istart_0=ifirst_0
          if (mod(istart_0+ic_1,2).ne.modit) istart_0=istart_0+1
          do ic_0=istart_0,ilast_0,2
            residual=
     &       matrix(ic_0,ic_1,-1,-1)*solution(ic_0-1,ic_1-1)
     &      +matrix(ic_0,ic_1, 0,-1)*solution(ic_0  ,ic_1-1)
     &      +matrix(ic_0,ic_1, 1,-1)*solution(ic_0+1,ic_1-1)
     &      +matrix(ic_0,ic_1,-1, 0)*solution(ic_0-1,ic_1  )
     &      +matrix(ic_0,ic_1, 0, 0)*solution(ic_0  ,ic_1  )
     &      +matrix(ic_0,ic_1, 1, 0)*solution(ic_0+1,ic_1  )
     &      +matrix(ic_0,ic_1,-1, 1)*solution(ic_0-1,ic_1+1)
     &      +matrix(ic_0,ic_1, 0, 1)*solution(ic_0  ,ic_1+1)
     &      +matrix(ic_0,ic_1, 1, 1)*solution(ic_0+1,ic_1+1)
     &      -rhs(ic_0,ic_1)
            residual_max=max(residual_max,abs(residual))
            solution(ic_0,ic_1)=solution(ic_0,ic_1)
     &                         -residual/matrix(ic_0,ic_1,0,0)
          enddo
        enddo
c       print *, "even solution"
c       call rbugcell(fi_0,fi_1,la_0,la_1,
c    &    ifirst_0,ifirst_1,ilast_0,ilast_1,solution)

        modit=1-modit
        do ic_1=ifirst_1,ilast_1
          istart_0=ifirst_0
          if (mod(istart_0+ic_1,2).ne.modit) istart_0=istart_0+1
          do ic_0=istart_0,ilast_0,2
            residual=
     &       matrix(ic_0,ic_1,-1,-1)*solution(ic_0-1,ic_1-1)
     &      +matrix(ic_0,ic_1, 0,-1)*solution(ic_0  ,ic_1-1)
     &      +matrix(ic_0,ic_1, 1,-1)*solution(ic_0+1,ic_1-1)
     &      +matrix(ic_0,ic_1,-1, 0)*solution(ic_0-1,ic_1  )
     &      +matrix(ic_0,ic_1, 0, 0)*solution(ic_0  ,ic_1  )
     &      +matrix(ic_0,ic_1, 1, 0)*solution(ic_0+1,ic_1  )
     &      +matrix(ic_0,ic_1,-1, 1)*solution(ic_0-1,ic_1+1)
     &      +matrix(ic_0,ic_1, 0, 1)*solution(ic_0  ,ic_1+1)
     &      +matrix(ic_0,ic_1, 1, 1)*solution(ic_0+1,ic_1+1)
     &      -rhs(ic_0,ic_1)
            residual_max=max(residual_max,abs(residual))
            solution(ic_0,ic_1)=solution(ic_0,ic_1)
     &                         -residual/matrix(ic_0,ic_1,0,0)
          enddo
        enddo
c       print *, "odd solution"
c       call rbugcell(fi_0,fi_1,la_0,la_1,
c    &    ifirst_0,ifirst_1,ilast_0,ilast_1,solution)
      else
c       print *, "odd iteration"
        modit=1
        do ic_1=ilast_1,ifirst_1,-1
          istart_0=ilast_0
          if (mod(istart_0+ic_1,2).ne.modit) istart_0=istart_0-1
          do ic_0=istart_0,ifirst_0+1,-2
            residual=
     &       matrix(ic_0,ic_1,-1,-1)*solution(ic_0-1,ic_1-1)
     &      +matrix(ic_0,ic_1, 0,-1)*solution(ic_0  ,ic_1-1)
     &      +matrix(ic_0,ic_1, 1,-1)*solution(ic_0+1,ic_1-1)
     &      +matrix(ic_0,ic_1,-1, 0)*solution(ic_0-1,ic_1  )
     &      +matrix(ic_0,ic_1, 0, 0)*solution(ic_0  ,ic_1  )
     &      +matrix(ic_0,ic_1, 1, 0)*solution(ic_0+1,ic_1  )
     &      +matrix(ic_0,ic_1,-1, 1)*solution(ic_0-1,ic_1+1)
     &      +matrix(ic_0,ic_1, 0, 1)*solution(ic_0  ,ic_1+1)
     &      +matrix(ic_0,ic_1, 1, 1)*solution(ic_0+1,ic_1+1)
     &      -rhs(ic_0,ic_1)
            residual_max=max(residual_max,abs(residual))
            solution(ic_0,ic_1)=solution(ic_0,ic_1)
     &                         -residual/matrix(ic_0,ic_1,0,0)
          enddo
        enddo
c       print *, "odd solution"
c       call rbugcell(fi_0,fi_1,la_0,la_1,
c    &    ifirst_0,ifirst_1,ilast_0,ilast_1,solution)

        modit=1-modit
        do ic_1=ilast_1,ifirst_1,-1
          istart_0=ilast_0
          if (mod(istart_0+ic_1,2).ne.modit) istart_0=istart_0-1
          do ic_0=istart_0,ifirst_0+1,-2
            residual=
     &       matrix(ic_0,ic_1,-1,-1)*solution(ic_0-1,ic_1-1)
     &      +matrix(ic_0,ic_1, 0,-1)*solution(ic_0  ,ic_1-1)
     &      +matrix(ic_0,ic_1, 1,-1)*solution(ic_0+1,ic_1-1)
     &      +matrix(ic_0,ic_1,-1, 0)*solution(ic_0-1,ic_1  )
     &      +matrix(ic_0,ic_1, 0, 0)*solution(ic_0  ,ic_1  )
     &      +matrix(ic_0,ic_1, 1, 0)*solution(ic_0+1,ic_1  )
     &      +matrix(ic_0,ic_1,-1, 1)*solution(ic_0-1,ic_1+1)
     &      +matrix(ic_0,ic_1, 0, 1)*solution(ic_0  ,ic_1+1)
     &      +matrix(ic_0,ic_1, 1, 1)*solution(ic_0+1,ic_1+1)
     &      -rhs(ic_0,ic_1)
            residual_max=max(residual_max,abs(residual))
            solution(ic_0,ic_1)=solution(ic_0,ic_1)
     &                         -residual/matrix(ic_0,ic_1,0,0)
          enddo
        enddo
c       print *, "even solution"
c       call rbugcell(fi_0,fi_1,la_0,la_1,
c    &    ifirst_0,ifirst_1,ilast_0,ilast_1,solution)
      endif
c     print *, "leaving gauss_seidel_red_black"
c     print *, "residual_max = ",residual_max
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

      subroutine sor(fi_0,fi_1,la_0,la_1,
     &  ifirst_0,ifirst_1,ilast_0,ilast_1, omega,residual_max,
     &  matrix,rhs,
     &  solution)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fi_0,fi_1,la_0,la_1,ifirst_0,ifirst_1,ilast_0,ilast_1
      double precision omega,residual_max
      double precision 
     &  matrix(fi_0:la_0,fi_1:la_1,-1:1,-1:1),
     &  rhs(fi_0:la_0,fi_1:la_1)
      double precision solution(fi_0:la_0,fi_1:la_1)

      integer ic_0,ic_1
      double precision residual
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      residual_max=0.d0
      do ic_1=ifirst_1,ilast_1
        do ic_0=ifirst_0,ilast_0
          residual=
     &       matrix(ic_0,ic_1,-1,-1)*solution(ic_0-1,ic_1-1)
     &      +matrix(ic_0,ic_1, 0,-1)*solution(ic_0  ,ic_1-1)
     &      +matrix(ic_0,ic_1, 1,-1)*solution(ic_0+1,ic_1-1)
     &      +matrix(ic_0,ic_1,-1, 0)*solution(ic_0-1,ic_1  )
     &      +matrix(ic_0,ic_1, 0, 0)*solution(ic_0  ,ic_1  )
     &      +matrix(ic_0,ic_1, 1, 0)*solution(ic_0+1,ic_1  )
     &      +matrix(ic_0,ic_1,-1, 1)*solution(ic_0-1,ic_1+1)
     &      +matrix(ic_0,ic_1, 0, 1)*solution(ic_0  ,ic_1+1)
     &      +matrix(ic_0,ic_1, 1, 1)*solution(ic_0+1,ic_1+1)
     &      -rhs(ic_0,ic_1)
          residual_max=max(residual_max,abs(residual))
          solution(ic_0,ic_1)=solution(ic_0,ic_1)
     &                       -omega*residual/matrix(ic_0,ic_1,0,0)
        enddo
      enddo
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

      subroutine iterative_improvement_dicf_setup(
     &  fi_0,fi_1,la_0,la_1,ifirst_0,ifirst_1,ilast_0,ilast_1,
     &  matrix,
     &  n,nnz,a,diag,col_ptr,row_ind)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fi_0,fi_1,la_0,la_1,ifirst_0,ifirst_1,ilast_0,ilast_1
      double precision matrix(fi_0:la_0,fi_1:la_1,-1:1,-1:1)
      integer n,nnz,col_ptr(1),row_ind(1)
      double precision diag(1),a(1)

      integer ic_0,ic_1,n_0,n_1
c     integer i,j
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "in iterative_improvement_dicf_setup"
c     do ic_1=ifirst_1,ilast_1
c       do ic_0=ifirst_0,ilast_0
c         print *, "matrix[",ic_0,",",ic_1,"] = ",
c    &      ((matrix(ic_0,ic_1,i),i=-1,1),j=-1,1)
c       enddo
c     enddo

      n_0=ilast_0-ifirst_0+1
      n_1=ilast_1-ifirst_1+1

      n=0
      nnz=0
      col_ptr(1)=1
      do ic_1=ifirst_1,ilast_1
        do ic_0=ifirst_0,ilast_0
          n=n+1
          diag(n)=matrix(ic_0,ic_1,0,0)
          if (ic_0.lt.ilast_0) then
            nnz=nnz+1
            row_ind(nnz)=n+1
            a(nnz)=matrix(ic_0,ic_1,1,0)
          endif
          if (ic_1.lt.ilast_1) then
            if (ic_0.gt.ifirst_0) then
              nnz=nnz+1
              row_ind(nnz)=n+n_0-1
              a(nnz)=matrix(ic_0,ic_1,-1,1)
            endif
            nnz=nnz+1
            row_ind(nnz)=n+n_0
            a(nnz)=matrix(ic_0,ic_1,0,1)
            if (ic_0.lt.ilast_0) then
              nnz=nnz+1
              row_ind(nnz)=n+n_0+1
              a(nnz)=matrix(ic_0,ic_1,1,1)
            endif
          endif
          col_ptr(n+1)=nnz+1
        enddo
      enddo

c     print *, "leaving iterative_improvement_dicf_setup"
c     print *, "n,nnz = ",n,nnz
c     do i=1,n
c       print *, "diag[",i,"] = ",diag(i)
c     enddo
c     do j=1,n+1
c       print *, "col_ptr[",j,"] = ",col_ptr(j)
c     enddo
c     do j=1,n
c       do i=col_ptr(j),col_ptr(j+1)-1
c         print *, "a[",row_ind(i),",",j,"] = ",a(i)
c       enddo
c     enddo
c     call flush(6)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

      subroutine iterative_improvement_dicf(fi_0,fi_1,la_0,la_1,
     &  ifirst_0,ifirst_1,ilast_0,ilast_1,n,nnz,residual_max, 
     &  a,diag,col_ptr,matrix,rhs,row_ind,
     &  solution,
     &  residual)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fi_0,fi_1,la_0,la_1,
     &  ifirst_0,ifirst_1,ilast_0,ilast_1,n,nnz
      double precision residual_max
      double precision a(nnz),diag(n), 
     &  matrix(fi_0:la_0,fi_1:la_1,-1:1,-1:1),
     &  rhs(fi_0:la_0,fi_1:la_1)
      integer col_ptr(n+1),row_ind(nnz)
      double precision solution(fi_0:la_0,fi_1:la_1)
      double precision residual(1)

      integer i,ic_0,ic_1,j,k
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "in iterative_improvement_dicf"
c     print *, "rhs"
c     call rbugcell(fi_0,fi_1,la_0,la_1,
c    &  ifirst_0,ifirst_1,ilast_0,ilast_1,rhs)
c     print *, "solution"
c     call rbugcell(fi_0,fi_1,la_0,la_1,
c    &  ifirst_0,ifirst_1,ilast_0,ilast_1,solution)
c     do j=-1,1
c       do i=-1,1
c         print *, "matrix(",i,",",j,")"
c         call rbugcell(fi_0,fi_1,la_0,la_1,
c    &      ifirst_0,ifirst_1,ilast_0,ilast_1,matrix(fi_0,fi_1,i,j))
c       enddo
c     enddo
      i=0
      residual_max=0.d0
      do ic_1=ifirst_1,ilast_1
        do ic_0=ifirst_0,ilast_0
          i=i+1
          residual(i)=
     &       matrix(ic_0,ic_1,-1,-1)*solution(ic_0-1,ic_1-1)
     &      +matrix(ic_0,ic_1, 0,-1)*solution(ic_0  ,ic_1-1)
     &      +matrix(ic_0,ic_1, 1,-1)*solution(ic_0+1,ic_1-1)
     &      +matrix(ic_0,ic_1,-1, 0)*solution(ic_0-1,ic_1  )
     &      +matrix(ic_0,ic_1, 0, 0)*solution(ic_0  ,ic_1  )
     &      +matrix(ic_0,ic_1, 1, 0)*solution(ic_0+1,ic_1  )
     &      +matrix(ic_0,ic_1,-1, 1)*solution(ic_0-1,ic_1+1)
     &      +matrix(ic_0,ic_1, 0, 1)*solution(ic_0  ,ic_1+1)
     &      +matrix(ic_0,ic_1, 1, 1)*solution(ic_0+1,ic_1+1)
     &      -rhs(ic_0,ic_1)
          residual_max=max(residual_max,abs(residual(i)))
        enddo
      enddo

c     forward-solve
      do j=1,n
        residual(j)=residual(j)/diag(j)
        do k=col_ptr(j),col_ptr(j+1)-1
          i=row_ind(k)
          residual(i)=residual(i)-a(k)*residual(j)
        enddo
      enddo
c     back-solve
      do j=n,1,-1
        do k=col_ptr(j),col_ptr(j+1)-1
          i=row_ind(k)
          residual(j)=residual(j)-a(k)*residual(i)
        enddo
        residual(j)=residual(j)/diag(j)
      enddo

      i=0
      do ic_1=ifirst_1,ilast_1
        do ic_0=ifirst_0,ilast_0
          i=i+1
          solution(ic_0,ic_1)=solution(ic_0,ic_1)-residual(i)
        enddo
      enddo
c     print *, "leaving iterative_improvement_dicf"
c     print *, "residual_max = ",residual_max
c     call flush(6)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

      subroutine matrix_multiply(fi_0,fi_1,la_0,la_1,
     &  ifirst_0,ifirst_1,ilast_0,ilast_1,
     &  a,x, ax)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fi_0,fi_1,la_0,la_1,
     &  ifirst_0,ifirst_1,ilast_0,ilast_1
      double precision 
     &  a(fi_0:la_0,fi_1:la_1,-1:1,-1:1),
     &  x(fi_0:la_0,fi_1:la_1)
      double precision 
     &  ax(fi_0:la_0,fi_1:la_1)

      integer ic_0,ic_1
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "in matrix_multiply"
c     print *, "x"
c     call rbugcell(fi_0,fi_1,la_0,la_1,
c    &  fi_0,fi_1,la_0,la_1,x)
c     do j=-1,1
c       do i=-1,1
c         print *, "a(",i,",",j,")"
c         call rbugcell(fi_0,fi_1,la_0,la_1,
c    &      ifirst_0,ifirst_1,ilast_0,ilast_1,a(fi_0,fi_1,i,j))
c       enddo
c     enddo
c     call flush(6)

      do ic_1=ifirst_1,ilast_1
        do ic_0=ifirst_0,ilast_0
          ax(ic_0,ic_1)=a(ic_0,ic_1,-1,-1)*x(ic_0-1,ic_1-1)
     &                 +a(ic_0,ic_1, 0,-1)*x(ic_0  ,ic_1-1)
     &                 +a(ic_0,ic_1, 1,-1)*x(ic_0+1,ic_1-1)
     &                 +a(ic_0,ic_1,-1, 0)*x(ic_0-1,ic_1  )
     &                 +a(ic_0,ic_1, 0, 0)*x(ic_0  ,ic_1  )
     &                 +a(ic_0,ic_1, 1, 0)*x(ic_0+1,ic_1  )
     &                 +a(ic_0,ic_1,-1, 1)*x(ic_0-1,ic_1+1)
     &                 +a(ic_0,ic_1, 0, 1)*x(ic_0  ,ic_1+1)
     &                 +a(ic_0,ic_1, 1, 1)*x(ic_0+1,ic_1+1)
        enddo
      enddo

c     print *, "in matrix_multiply"
c     print *, "ax"
c     call rbugcell(fi_0,fi_1,la_0,la_1,
c    &  ifirst_0,ifirst_1,ilast_0,ilast_1,ax)
c     call flush(6)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

      subroutine jacobi_preconditioner(fi_0,fi_1,la_0,la_1,
     &  ifirst_0,ifirst_1,ilast_0,ilast_1,
     &  matrix,rhs,
     &  solution)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fi_0,fi_1,la_0,la_1,ifirst_0,ifirst_1,ilast_0,ilast_1
      double precision 
     &  matrix(fi_0:la_0,fi_1:la_1,-1:1,-1:1),
     &  rhs(fi_0:la_0,fi_1:la_1)
      double precision solution(fi_0:la_0,fi_1:la_1)

      integer ic_0,ic_1
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      do ic_1=ifirst_1,ilast_1
        do ic_0=ifirst_0,ilast_0
          solution(ic_0,ic_1)=rhs(ic_0,ic_1)/matrix(ic_0,ic_1,0,0)
        enddo
      enddo
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

      subroutine block_jacobi_preconditioner(fi_0,fi_1,la_0,la_1,
     &  ifirst_0,ifirst_1,ilast_0,ilast_1,
     &  matrix,rhs,
     &  diagonal_block,solution)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fi_0,fi_1,la_0,la_1,ifirst_0,ifirst_1,ilast_0,ilast_1
      double precision 
     &  matrix(fi_0:la_0,fi_1:la_1,-1:1,-1:1),
     &  rhs(fi_0:la_0,fi_1:la_1)
      double precision
     &  diagonal_block(fi_0:la_0,0:1),
     &  solution(fi_0:la_0,fi_1:la_1)

      integer ic_0,ic_1
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      do ic_1=ifirst_1,ilast_1
        do ic_0=ifirst_0,ilast_0
          diagonal_block(ic_0, 0)=matrix(ic_0,ic_1,0,0)
          diagonal_block(ic_0, 1)=matrix(ic_0,ic_1,1,0)
          solution(ic_0,ic_1)=rhs(ic_0,ic_1)
        enddo
        do ic_0=ifirst_0,ilast_0-1
          diagonal_block(ic_0,0)=1.d0/sqrt(diagonal_block(ic_0,0))
          diagonal_block(ic_0,1)=diagonal_block(ic_0,1)
     &                          *diagonal_block(ic_0,0)
          diagonal_block(ic_0+1,0)=diagonal_block(ic_0+1,0)
     &                            -diagonal_block(ic_0,1)**2
          solution(ic_0,ic_1)=solution(ic_0,ic_1)*diagonal_block(ic_0,0)
          solution(ic_0+1,ic_1)=solution(ic_0+1,ic_1)
     &                       -diagonal_block(ic_0,1)*solution(ic_0,ic_1)
        enddo
        solution(ilast_0,ic_1)=solution(ilast_0,ic_1)
     &                        *diagonal_block(ilast_0,0)
        do ic_0=ilast_0-1,ifirst_0,-1
          solution(ic_0,ic_1)=(solution(ic_0,ic_1)
     &      -diagonal_block(ic_0,1)*solution(ic_0+1,ic_1))
     &      *diagonal_block(ic_0,0)
        enddo
      enddo
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

      subroutine dicf_preconditioner(fi_0,fi_1,la_0,la_1,
     &  ifirst_0,ifirst_1,ilast_0,ilast_1,n,nnz, 
     &  a,diag,col_ptr,rhs,row_ind,
     &  solution,soln)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fi_0,fi_1,la_0,la_1,
     &  ifirst_0,ifirst_1,ilast_0,ilast_1,n,nnz
      double precision a(nnz),diag(n), 
     &  rhs(fi_0:la_0,fi_1:la_1)
      integer col_ptr(n+1),row_ind(nnz)
      double precision solution(fi_0:la_0,fi_1:la_1),soln(n)

      integer i,ic_0,ic_1,j,k
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      i=0
      do ic_1=ifirst_1,ilast_1
        do ic_0=ifirst_0,ilast_0
          i=i+1
          soln(i)=rhs(ic_0,ic_1)
        enddo
      enddo

c     forward-solve
      do j=1,n
        soln(j)=soln(j)/diag(j)
        do k=col_ptr(j),col_ptr(j+1)-1
          i=row_ind(k)
          soln(i)=soln(i)-a(k)*soln(j)
        enddo
      enddo
c     back-solve
      do j=n,1,-1
        do k=col_ptr(j),col_ptr(j+1)-1
          i=row_ind(k)
          soln(j)=soln(j)-a(k)*soln(i)
        enddo
        soln(j)=soln(j)/diag(j)
      enddo

      i=0
      do ic_1=ifirst_1,ilast_1
        do ic_0=ifirst_0,ilast_0
          i=i+1
          solution(ic_0,ic_1)=soln(i)
        enddo
      enddo
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end 

      subroutine vector_setup(fc_0,fc_1,lc_0,lc_1,
     &  start_0,start_1,vector)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fc_0,fc_1,lc_0,lc_1,start_0,start_1
      double precision vector(fc_0:lc_0,fc_1:lc_1)
      integer ic_0,ic_1
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "in vector_setup"
c     print *, "fc = ",fc_0,fc_1
c     print *, "lc = ",lc_0,lc_1
c     print *, "start = ",start_0,start_1
c     call flush(6)

      do ic_1=fc_1+start_1,lc_1,3
        do ic_0=fc_0+start_0,lc_0,3
          vector(ic_0,ic_1)=1.d0
        enddo
      enddo

c     print *, "leaving vector_setup"
c     call rbugcell(fc_0,fc_1,lc_0,lc_1,
c    &  fc_0,fc_1,lc_0,lc_1,vector)
c     call flush(6)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

      subroutine matrix_setup(fc_0,fc_1,lc_0,lc_1,
     &  start_0,start_1,vector,matrix)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fc_0,fc_1,lc_0,lc_1,start_0,start_1
      double precision vector(fc_0:lc_0,fc_1:lc_1)
      double precision matrix(fc_0:lc_0,fc_1:lc_1,-1:1,-1:1)
      integer ic_0,ic_1
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     assumes that boundary values are defined
      do ic_1=fc_1+start_1,lc_1,3
        do ic_0=fc_0+start_0,lc_0,3
          matrix(ic_0,ic_1,-1,-1)=vector(ic_0-1,ic_1-1)
          matrix(ic_0,ic_1, 0,-1)=vector(ic_0  ,ic_1-1)
          matrix(ic_0,ic_1, 1,-1)=vector(ic_0+1,ic_1-1)
          matrix(ic_0,ic_1,-1, 0)=vector(ic_0-1,ic_1  )
          matrix(ic_0,ic_1, 0, 0)=vector(ic_0  ,ic_1  )
          matrix(ic_0,ic_1, 1, 0)=vector(ic_0+1,ic_1  )
          matrix(ic_0,ic_1,-1, 1)=vector(ic_0-1,ic_1+1)
          matrix(ic_0,ic_1, 0, 1)=vector(ic_0  ,ic_1+1)
          matrix(ic_0,ic_1, 1, 1)=vector(ic_0+1,ic_1+1)
        enddo
      enddo
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

      subroutine prolong_fem(fc_0,fc_1,lc_0,lc_1,
     &  ff_0,ff_1,lf_0,lf_1,
     &  vectorc,vectorf)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fc_0,fc_1,lc_0,lc_1,
     &  ff_0,ff_1,lf_0,lf_1
      double precision 
     &  vectorc(fc_0:lc_0,fc_1:lc_1)
      double precision 
     &  vectorf(ff_0:lf_0,ff_1:lf_1)
      integer ic_0,ic_1,if_0,if_1
      double precision sum
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     assumes Dirichlet boundary values
c     Briggs, p. 39:
      do if_1=ff_1+2,lf_1-1,2
        ic_1=if_1/2
        do if_0=ff_0+2,lf_0-1,2
          ic_0=if_0/2
          vectorf(if_0,if_1)=vectorc(ic_0,ic_1)
        enddo
        do if_0=ff_0+1,lf_0-1,2
          ic_0=(if_0-1)/2
          sum=0.d0
          if (ic_0.gt.fc_0) sum=sum+vectorc(ic_0,ic_1)
          if (ic_0+1.lt.lc_0) sum=sum+vectorc(ic_0+1,ic_1)
          vectorf(if_0,if_1)=0.5d0*sum
        enddo
      enddo
      do if_1=ff_1+1,lf_1-1,2
        ic_1=(if_1-1)/2
        do if_0=ff_0+2,lf_0-1,2
          ic_0=if_0/2
          sum=0.d0
          if (ic_1.gt.fc_1) sum=sum+vectorc(ic_0,ic_1)
          if (ic_1+1.lt.lc_1) sum=sum+vectorc(ic_0,ic_1+1)
          vectorf(if_0,if_1)=0.5d0*sum
        enddo
        do if_0=ff_0+1,lf_0-1,2
          ic_0=(if_0-1)/2
          sum=0.d0
          if (ic_0.gt.fc_0) then
            if (ic_1.gt.fc_1) sum=sum+vectorc(ic_0,ic_1)
            if (ic_1+1.lt.lc_1) sum=sum+vectorc(ic_0,ic_1+1)
          endif
          if (ic_0+1.lt.lc_0) then
            if (ic_1.gt.fc_1) sum=sum+vectorc(ic_0+1,ic_1)
            if (ic_1+1.lt.lc_1) sum=sum+vectorc(ic_0+1,ic_1+1)
          endif
          vectorf(if_0,if_1)=0.25d0*sum
        enddo
      enddo
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

      subroutine restrict_fem(fc_0,fc_1,lc_0,lc_1,
     &  ff_0,ff_1,lf_0,lf_1,
     &  vectorf,vectorc)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fc_0,fc_1,lc_0,lc_1,
     &  ff_0,ff_1,lf_0,lf_1
      double precision 
     &  vectorf(ff_0:lf_0,ff_1:lf_1)
      double precision 
     &  vectorc(fc_0:lc_0,fc_1:lc_1)
      integer ic_0,ic_1,if_0,if_1
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     assumes Dirichlet boundary values
c     Briggs, p. 40:
      do ic_1=fc_1+1,lc_1-1
        if_1=2*ic_1
        do ic_0=fc_0+1,lc_0-1
          if_0=2*ic_0
c         vectorc(ic_0,ic_1)=0.25d0*(vectorf(if_0,if_1)
c    &      +0.5d0*(vectorf(if_0-1,if_1  )+vectorf(if_0+1,if_1  )
c    &             +vectorf(if_0  ,if_1-1)+vectorf(if_0  ,if_1+1))
c    &      +0.25d0*(vectorf(if_0-1,if_1-1)+vectorf(if_0+1,if_1-1)
c    &              +vectorf(if_0-1,if_1+1)+vectorf(if_0+1,if_1+1)))
          vectorc(ic_0,ic_1)=vectorf(if_0,if_1)
     &      +0.5d0*(vectorf(if_0-1,if_1  )+vectorf(if_0+1,if_1  )
     &             +vectorf(if_0  ,if_1-1)+vectorf(if_0  ,if_1+1))
     &      +0.25d0*(vectorf(if_0-1,if_1-1)+vectorf(if_0+1,if_1-1)
     &              +vectorf(if_0-1,if_1+1)+vectorf(if_0+1,if_1+1))
        enddo
      enddo
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

      subroutine setup_am(fc_0,fc_1,lc_0,lc_1,
     &  ff_0,ff_1,lf_0,lf_1,
     & matf, prolongcell,prolongside0,prolongside1)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fc_0,fc_1,lc_0,lc_1, ! cell indices
     &  ff_0,ff_1,lf_0,lf_1
      double precision 
     &  matf(ff_0:lf_0+1,ff_1:lf_1+1,-1:1,-1:1) ! corner-centered
      double precision 
     &  prolongcell(fc_0:lc_0,fc_1:lc_1,2,2),
     &  prolongside0(fc_0:lc_0+1,fc_1:lc_1,2),
     &  prolongside1(fc_1:lc_1+1,fc_0:lc_0,2)
      integer ic_0,ic_1,if_0,if_1
      double precision den,sum
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "in setup_am"
c     print *, "fc = ",fc_0,fc_1
c     print *, "lc = ",lc_0,lc_1
c     print *, "ff = ",ff_0,ff_1
c     print *, "lf = ",lf_0,lf_1
c     do j=-1,1
c       do i=-1,1
c         print *, "matf(",i,",",j,") = "
c         call rbugcorn(ff_0,ff_1,lf_0,lf_1,
c    &      ff_0,ff_1,lf_0,lf_1,matf(ff_0,ff_1,i,j))
c       enddo
c     enddo
c     call flush(6)

c     print *, " "
c     call flush(6)
      do ic_1=fc_1,lc_1 ! coarse cell index
        if_1=2*ic_1+1 ! fine corner index
c       print *, "ic_1,if_1 = ",ic_1,if_1
c       call flush(6)

        ic_0=fc_0
c       print *, "   ic_0 = ",ic_0
c       call flush(6)
        prolongside0(ic_0,ic_1,1)=0.d0
        prolongside0(ic_0,ic_1,2)=0.d0
        if (ic_1.le.fc_1) then
          do ic_0=fc_0+1,lc_0 ! coarse side index
            if_0=2*ic_0 ! fine corner index
c           print *, "   ic_0,if_0 = ",ic_0,if_0
c           call flush(6)
            sum=
     &      matf(if_0,if_1,-1,1)+matf(if_0,if_1,0,1)+matf(if_0,if_1,1,1)
     &     +matf(if_0,if_1,-1,0)                    +matf(if_0,if_1,1,0)
            prolongside0(ic_0,ic_1,1)=0.d0
            prolongside0(ic_0,ic_1,2)=-sum/matf(if_0,if_1,0,0)
          enddo
        else if (ic_1.lt.lc_1) then
          do ic_0=fc_0+1,lc_0 ! coarse side index
            if_0=2*ic_0 ! fine corner index
c           print *, "   ic_0,if_0 = ",ic_0,if_0
c           call flush(6)
            den=matf(if_0,if_1,0,-1)+matf(if_0,if_1,0,1)
            sum=
     &   matf(if_0,if_1,-1, 1)+matf(if_0,if_1,0, 1)+matf(if_0,if_1,1, 1)
     &  +matf(if_0,if_1,-1, 0)                     +matf(if_0,if_1,1, 1)
     &  +matf(if_0,if_1,-1,-1)+matf(if_0,if_1,0,-1)+matf(if_0,if_1,1,-1)
            sum=sum/matf(if_0,if_1,0,0)
            prolongside0(ic_0,ic_1,1)=-sum*matf(if_0,if_1,0,-1)/den
            prolongside0(ic_0,ic_1,2)=-sum*matf(if_0,if_1,0, 1)/den
          enddo
        else
          do ic_0=fc_0+1,lc_0 ! coarse side index
            if_0=2*ic_0 ! fine corner index
c           print *, "   ic_0,if_0 = ",ic_0,if_0
c           call flush(6)
            sum=
     &   matf(if_0,if_1,-1, 0)                     +matf(if_0,if_1,1, 1)
     &  +matf(if_0,if_1,-1,-1)+matf(if_0,if_1,0,-1)+matf(if_0,if_1,1,-1)
            prolongside0(ic_0,ic_1,1)=-sum/matf(if_0,if_1,0,0)
            prolongside0(ic_0,ic_1,2)=0.d0
          enddo
        endif
        ic_0=lc_0+1
c       print *, "   ic_0 = ",ic_0
c       call flush(6)
        prolongside0(ic_0,ic_1,1)=0.d0
        prolongside0(ic_0,ic_1,2)=0.d0
      enddo

c     print *, " "
c     call flush(6)
      do ic_0=fc_0,lc_0 ! coarse cell index
        if_0=2*ic_0+1 ! fine corner index
c       print *, "ic_0,if_0 = ",ic_0,if_0
c       call flush(6)

        ic_1=fc_1
c       print *, "   ic_1 = ",ic_1
c       call flush(6)
        prolongside1(ic_1,ic_0,1)=0.d0
        prolongside1(ic_1,ic_0,2)=0.d0
        if (ic_0.le.fc_0) then
          do ic_1=fc_1+1,lc_1
            if_1=2*ic_1
c           print *, "   ic_1,if_1 = ",ic_1,if_1
c           call flush(6)
            sum=
     &                         matf(if_0,if_1,0, 1)+matf(if_0,if_1,1, 1)
     &                                             +matf(if_0,if_1,1, 0)
     &                        +matf(if_0,if_1,0,-1)+matf(if_0,if_1,1,-1)
            prolongside1(ic_1,ic_0,1)=0.d0
            prolongside1(ic_1,ic_0,2)=-sum/matf(if_0,if_1,0,0)
          enddo
        else if (ic_0.lt.lc_0) then
          do ic_1=fc_1+1,lc_1
            if_1=2*ic_1
c           print *, "   ic_1,if_1 = ",ic_1,if_1
c           call flush(6)
            den=matf(if_0,if_1,-1,0)+matf(if_0,if_1,1,0)
            sum=
     &   matf(if_0,if_1,-1, 1)+matf(if_0,if_1,0, 1)+matf(if_0,if_1,1, 1)
     &  +matf(if_0,if_1,-1, 0)                     +matf(if_0,if_1,1, 1)
     &  +matf(if_0,if_1,-1,-1)+matf(if_0,if_1,0,-1)+matf(if_0,if_1,1,-1)
            sum=sum/matf(if_0,if_1,0,0)
            prolongside1(ic_1,ic_0,1)=-sum*matf(if_0,if_1,-1,0)/den
            prolongside1(ic_1,ic_0,2)=-sum*matf(if_0,if_1, 1,0)/den
          enddo
        else
          do ic_1=fc_1+1,lc_1
            if_1=2*ic_1
c           print *, "   ic_1,if_1 = ",ic_1,if_1
c           call flush(6)
            sum=
     &      matf(if_0,if_1,-1, 1)+matf(if_0,if_1,0, 1)
     &     +matf(if_0,if_1,-1, 0)                    
     &     +matf(if_0,if_1,-1,-1)+matf(if_0,if_1,0,-1)
            prolongside1(ic_1,ic_0,1)=-sum/matf(if_0,if_1,0,0)
            prolongside1(ic_1,ic_0,2)=0.d0
          enddo
        endif
        ic_1=lc_1+1
c       print *, "   ic_1 = ",ic_1
c       call flush(6)
        prolongside1(ic_1,ic_0,1)=0.d0
        prolongside1(ic_1,ic_0,2)=0.d0
      enddo

      do ic_1=fc_1,lc_1 ! coarse cell index
        if_1=2*ic_1+1 ! fine corner index
c       print *, "ic_1,if_1 = ",ic_1,if_1
c       call flush(6)
        if (ic_1.le.fc_1) then
          do ic_0=fc_0,lc_0 ! coarse cell index
            prolongcell(ic_0,ic_1,1,1)=0.d0
            prolongcell(ic_0,ic_1,2,1)=0.d0
            if_0=2*ic_0+1 ! fine corner index
c           print *, "ic_0,if_0 = ",ic_0,if_0
c           call flush(6)
            if (ic_0.le.fc_0) then
              sum=
     &                        matf(if_0,if_1,0, 1)+matf(if_0,if_1,1, 1)
     &                                            +matf(if_0,if_1,1, 1)
              sum=sum/matf(if_0,if_1,0,0)
              den=
     &                       matf(if_0,if_1,0, 1)
     &                                           +matf(if_0,if_1,1, 1)
              prolongcell(ic_0,ic_1,1,2)=0.d0
              prolongcell(ic_0,ic_1,2,2)=
     &          -sum*(matf(if_0,if_1, 1, 0)*prolongside0(ic_0+1,ic_1,2)
     &               +matf(if_0,if_1, 0, 1)*prolongside1(ic_1+1,ic_0,2))
     &          /den
            else if (ic_0.lt.lc_0) then
              sum=
     &  matf(if_0,if_1,-1, 1)+matf(if_0,if_1,0, 1)+matf(if_0,if_1,1, 1)
     & +matf(if_0,if_1,-1, 0)                     +matf(if_0,if_1,1, 1)
              sum=sum/matf(if_0,if_1,0,0)
              den=
     &                        matf(if_0,if_1,0, 1)
     & +matf(if_0,if_1,-1, 0)                     +matf(if_0,if_1,1, 1)
              prolongcell(ic_0,ic_1,1,2)=
     &          -sum*(matf(if_0,if_1,-1, 0)*prolongside0(ic_0  ,ic_1,2)
     &               +matf(if_0,if_1, 0, 1)*prolongside1(ic_1+1,ic_0,1))
     &          /den
              prolongcell(ic_0,ic_1,2,2)=
     &          -sum*(matf(if_0,if_1, 1, 0)*prolongside0(ic_0+1,ic_1,2)
     &               +matf(if_0,if_1, 0, 1)*prolongside1(ic_1+1,ic_0,2))
     &          /den
            else
              sum=
     &    matf(if_0,if_1,-1, 1)+matf(if_0,if_1,0, 1)
     &   +matf(if_0,if_1,-1, 0)
              sum=sum/matf(if_0,if_1,0,0)
              den=
     &                          matf(if_0,if_1,0, 1)
     &   +matf(if_0,if_1,-1, 0)
              prolongcell(ic_0,ic_1,1,2)=
     &          -sum*(matf(if_0,if_1,-1, 0)*prolongside0(ic_0  ,ic_1,2)
     &               +matf(if_0,if_1, 0, 1)*prolongside1(ic_1+1,ic_0,1))
     &          /den
              prolongcell(ic_0,ic_1,2,2)=0.d0
            endif
          enddo
        else if (ic_1.lt.lc_1) then
          do ic_0=fc_0,lc_0 ! coarse cell index
            if_0=2*ic_0+1 ! fine corner index
c           print *, "ic_0,if_0 = ",ic_0,if_0
c           call flush(6)
            if (ic_0.le.fc_0) then
              sum=
     &                     matf(if_0,if_1,0, 1)+matf(if_0,if_1,1, 1)
     &                                         +matf(if_0,if_1,1, 1)
     &                    +matf(if_0,if_1,0,-1)+matf(if_0,if_1,1,-1)
              sum=sum/matf(if_0,if_1,0,0)
              den=
     &                     matf(if_0,if_1,0, 1)
     &                                         +matf(if_0,if_1,1, 1)
     &                    +matf(if_0,if_1,0,-1)
              prolongcell(ic_0,ic_1,1,1)=0.d0
              prolongcell(ic_0,ic_1,2,1)=
     &          -sum*(matf(if_0,if_1, 1, 0)*prolongside0(ic_0+1,ic_1,1)
     &               +matf(if_0,if_1, 0,-1)*prolongside1(ic_1  ,ic_0,2))
     &          /den
              prolongcell(ic_0,ic_1,1,2)=0.d0
              prolongcell(ic_0,ic_1,2,2)=
     &          -sum*(matf(if_0,if_1, 1, 0)*prolongside0(ic_0+1,ic_1,2)
     &               +matf(if_0,if_1, 0, 1)*prolongside1(ic_1+1,ic_0,2))
     &          /den
            else if (ic_0.lt.lc_0) then
              sum=
     &  matf(if_0,if_1,-1, 1)+matf(if_0,if_1,0, 1)+matf(if_0,if_1,1, 1)
     & +matf(if_0,if_1,-1, 0)                     +matf(if_0,if_1,1, 1)
     & +matf(if_0,if_1,-1,-1)+matf(if_0,if_1,0,-1)+matf(if_0,if_1,1,-1)
              sum=sum/matf(if_0,if_1,0,0)
              den=
     &                        matf(if_0,if_1,0, 1)
     & +matf(if_0,if_1,-1, 0)                     +matf(if_0,if_1,1, 1)
     &                       +matf(if_0,if_1,0,-1)
              prolongcell(ic_0,ic_1,1,1)=
     &          -sum*(matf(if_0,if_1,-1, 0)*prolongside0(ic_0  ,ic_1,1)
     &               +matf(if_0,if_1, 0,-1)*prolongside1(ic_1  ,ic_0,1))
     &          /den
              prolongcell(ic_0,ic_1,2,1)=
     &          -sum*(matf(if_0,if_1, 1, 0)*prolongside0(ic_0+1,ic_1,1)
     &               +matf(if_0,if_1, 0,-1)*prolongside1(ic_1  ,ic_0,2))
     &          /den
              prolongcell(ic_0,ic_1,1,2)=
     &          -sum*(matf(if_0,if_1,-1, 0)*prolongside0(ic_0  ,ic_1,2)
     &               +matf(if_0,if_1, 0, 1)*prolongside1(ic_1+1,ic_0,1))
     &          /den
              prolongcell(ic_0,ic_1,2,2)=
     &          -sum*(matf(if_0,if_1, 1, 0)*prolongside0(ic_0+1,ic_1,2)
     &               +matf(if_0,if_1, 0, 1)*prolongside1(ic_1+1,ic_0,2))
     &          /den
            else
              sum=
     &  matf(if_0,if_1,-1, 1)+matf(if_0,if_1,0, 1)
     & +matf(if_0,if_1,-1, 0)
     & +matf(if_0,if_1,-1,-1)+matf(if_0,if_1,0,-1)
              sum=sum/matf(if_0,if_1,0,0)
              den=
     &                        matf(if_0,if_1,0, 1)
     & +matf(if_0,if_1,-1, 0)
     &                       +matf(if_0,if_1,0,-1)
              prolongcell(ic_0,ic_1,1,1)=
     &          -sum*(matf(if_0,if_1,-1, 0)*prolongside0(ic_0  ,ic_1,1)
     &               +matf(if_0,if_1, 0,-1)*prolongside1(ic_1  ,ic_0,1))
     &          /den
              prolongcell(ic_0,ic_1,2,1)=0.d0
              prolongcell(ic_0,ic_1,1,2)=
     &          -sum*(matf(if_0,if_1,-1, 0)*prolongside0(ic_0  ,ic_1,2)
     &               +matf(if_0,if_1, 0, 1)*prolongside1(ic_1+1,ic_0,1))
     &          /den
              prolongcell(ic_0,ic_1,2,2)=0.d0
            endif
          enddo
        else
          do ic_0=fc_0,lc_0 ! coarse cell index
            if_0=2*ic_0+1 ! fine corner index
c           print *, "ic_0,if_0 = ",ic_0,if_0
c           call flush(6)
            prolongcell(ic_0,ic_1,1,2)=0.d0
            prolongcell(ic_0,ic_1,2,2)=0.d0
            if (ic_0.le.fc_0) then
              sum=
     &                                             matf(if_0,if_1,1, 1)
     &                       +matf(if_0,if_1,0,-1)+matf(if_0,if_1,1,-1)
              sum=sum/matf(if_0,if_1,0,0)
              den=
     &                                             matf(if_0,if_1,1, 1)
     &                       +matf(if_0,if_1,0,-1)
              prolongcell(ic_0,ic_1,1,1)=0.d0
              prolongcell(ic_0,ic_1,2,1)=
     &          -sum*(matf(if_0,if_1, 1, 0)*prolongside0(ic_0+1,ic_1,1)
     &               +matf(if_0,if_1, 0,-1)*prolongside1(ic_1  ,ic_0,2))
     &          /den
            else if (ic_0.lt.lc_0) then
              sum=
     &  matf(if_0,if_1,-1, 0)                     +matf(if_0,if_1,1, 1)
     & +matf(if_0,if_1,-1,-1)+matf(if_0,if_1,0,-1)+matf(if_0,if_1,1,-1)
              sum=sum/matf(if_0,if_1,0,0)
              den=
     &  matf(if_0,if_1,-1, 0)                     +matf(if_0,if_1,1, 1)
     &                       +matf(if_0,if_1,0,-1)
              prolongcell(ic_0,ic_1,1,1)=
     &          -sum*(matf(if_0,if_1,-1, 0)*prolongside0(ic_0  ,ic_1,1)
     &               +matf(if_0,if_1, 0,-1)*prolongside1(ic_1  ,ic_0,1))
     &          /den
              prolongcell(ic_0,ic_1,2,1)=
     &          -sum*(matf(if_0,if_1, 1, 0)*prolongside0(ic_0+1,ic_1,1)
     &               +matf(if_0,if_1, 0,-1)*prolongside1(ic_1  ,ic_0,2))
     &          /den
            else
              sum=
     &  matf(if_0,if_1,-1, 0)
     & +matf(if_0,if_1,-1,-1)+matf(if_0,if_1,0,-1)
              sum=sum/matf(if_0,if_1,0,0)
              den=
     &  matf(if_0,if_1,-1, 0)
     &                       +matf(if_0,if_1,0,-1)
              prolongcell(ic_0,ic_1,1,1)=
     &          -sum*(matf(if_0,if_1,-1, 0)*prolongside0(ic_0  ,ic_1,1)
     &               +matf(if_0,if_1, 0,-1)*prolongside1(ic_1  ,ic_0,1))
     &          /den
              prolongcell(ic_0,ic_1,2,1)=0.d0
            endif
          enddo
        endif
      enddo

c     print *, "leaving setup_am"
c     do i=1,2
c       print *, "prolongside0(*,",i,") = "
c       call rbugside0(fc_0,fc_1,lc_0,lc_1,
c    &    fc_0,fc_1,lc_0,lc_1,prolongside0(fc_0,fc_1,i))
c     enddo
c     do i=1,2
c       print *, "prolongside1(*,",i,") = "
c       call rbugside1(fc_0,fc_1,lc_0,lc_1,
c    &    fc_0,fc_1,lc_0,lc_1,prolongside1(fc_1,fc_0,i))
c     enddo
c     do j=1,2
c       do i=1,2
c         print *, "prolongcell(*,",i,",",j,") = "
c         call rbugcell(fc_0,fc_1,lc_0,lc_1,
c    &      fc_0,fc_1,lc_0,lc_1,prolongcell(fc_0,fc_1,i,j))
c       enddo
c     enddo
c     call flush(6)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

      subroutine prolong_am(fc_0,fc_1,lc_0,lc_1,
     &  ff_0,ff_1,lf_0,lf_1,
     &  prolongcell,prolongside0,prolongside1,vectorc,vectorf)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c***********************************************************************
c  Copyright 2006 John A. Trangenstein
c
c  This software is made available for research and instructional use 
c  only. 
c  You may copy and use this software without charge for these 
c  non-commercial purposes, provided that the copyright notice and 
c  associated text is reproduced on all copies.  
c  For all other uses (including distribution of modified versions), 
c  please contact the author at
c    John A. Trangenstein
c    Department of Mathematics
c    Duke University
c    Durham, NC 27708-0320
c    USA
c  or
c    johnt@math.duke.edu
c  
c  This software is made available "as is" without any assurance that it
c  is completely correct, or that it will work for your purposes.  
c  Use the software at your own risk.
c***********************************************************************
      double precision roundoff,small,huge,undefind
      common/machine/roundoff,small,huge,undefind

      double precision e,pi,rt2
      double precision half,one,two,zero
      integer debug_on
      common/const/e,pi,rt2, half,one,two,zero,debug_on
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fc_0,fc_1,lc_0,lc_1,
     &  ff_0,ff_1,lf_0,lf_1
      double precision 
     &  prolongcell(fc_0:lc_0,fc_1:lc_1,2,2),
     &  prolongside0(fc_0:lc_0+1,fc_1:lc_1,2),
     &  prolongside1(fc_1:lc_1+1,fc_0:lc_0,2),
     &  vectorc(fc_0:lc_0+1,fc_1:lc_1+1)
      double precision 
     &  vectorf(ff_0:lf_0+1,ff_1:lf_1+1)
      integer ic_0,ic_1,if_0,if_1
c     integer i,j
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     if (debug_on.eq.1) then
c       print *, "in prolong_am"
c       print *, "fc = ",fc_0,fc_1
c       print *, "lc = ",lc_0,lc_1
c       print *, "ff = ",ff_0,ff_1
c       print *, "lf = ",lf_0,lf_1
c       print *, "vectorc = "
c       call rbugcorn(fc_0,fc_1,lc_0,lc_1,
c    &    fc_0,fc_1,lc_0,lc_1,vectorc)
c       do j=1,2
c         do i=1,2
c           print *, "prolongcell(",i,",",j,") = "
c           call rbugcell(fc_0,fc_1,lc_0,lc_1,
c    &        fc_0,fc_1,lc_0,lc_1,prolongcell(fc_0,fc_1,i,j))
c         enddo
c       enddo
c       do i=1,2
c         print *, "prolongside0(*,",i,") = "
c         call rbugside0(fc_0,fc_1,lc_0,lc_1,
c    &      fc_0,fc_1,lc_0,lc_1,prolongside0(fc_0,fc_1,i))
c       enddo
c       do i=1,2
c         print *, "prolongside1(*,",i,") = "
c         call rbugside1(fc_0,fc_1,lc_0,lc_1,
c    &      fc_0,fc_1,lc_0,lc_1,prolongside1(fc_1,fc_0,i))
c       enddo
c       call flush(6)
c     endif

c     if (debug_on.eq.1) then
c       print *, "even,even"
c       call flush(6)
c     endif
c     assumes Dirichlet boundary values
      do ic_1=fc_1+1,lc_1 ! corner index
        if_1=ic_1*2
        do ic_0=fc_0+1,lc_0 ! corner index
          if_0=ic_0*2
          vectorf(if_0,if_1)=vectorc(ic_0,ic_1)
c         if (debug_on.eq.1) then
c           print *, "vectorf[",if_0,",",if_1,"] = ",vectorf(if_0,if_1)
c           call flush(6)
c         endif
        enddo
      enddo

      if (debug_on.eq.1) then
        print *, "odd,even"
        call flush(6)
      endif
      do ic_1=fc_1+1,lc_1 ! corner index
        if_1=ic_1*2       ! corner index
        do ic_0=fc_0,lc_0 ! cell index
          if_0=ic_0*2+1   ! corner index
          vectorf(if_0,if_1)= ! center of coarse side ic_1,ic_0
     &      prolongside1(ic_1,ic_0,1)*vectorf(if_0-1,if_1)
     &     +prolongside1(ic_1,ic_0,2)*vectorf(if_0+1,if_1)
c         if (debug_on.eq.1) then
c           print *, "vectorf[",if_0,",",if_1,"] = ",vectorf(if_0,if_1)
c           call flush(6)
c         endif
        enddo
      enddo

      if (debug_on.eq.1) then
        print *, "even,odd"
        call flush(6)
      endif
      do ic_0=fc_0+1,lc_0 ! corner index
        if_0=ic_0*2       ! corner index
        do ic_1=fc_1,lc_1 ! cell index
          if_1=ic_1*2+1   ! corner index
          vectorf(if_0,if_1)= ! center of coarse side ic_0,ic_1
     &      prolongside0(ic_0,ic_1,1)*vectorf(if_0,if_1-1)
     &     +prolongside0(ic_0,ic_1,2)*vectorf(if_0,if_1+1)
c         if (debug_on.eq.1) then
c           print *, "vectorf[",if_0,",",if_1,"] = ",vectorf(if_0,if_1)
c           call flush(6)
c         endif
        enddo
      enddo

      if (debug_on.eq.1) then
        print *, "odd,odd"
        call flush(6)
      endif
      do ic_1=fc_1,lc_1   ! cell index
        if_1=ic_1*2+1     ! corner index
        do ic_0=fc_0,lc_0 ! cell index
          if_0=ic_0*2+1   ! corner index
          vectorf(if_0,if_1)= ! center of coarse cell ic_0,ic_1
     &      prolongcell(ic_0,ic_1,1,1)*vectorf(if_0-1,if_1-1)
     &     +prolongcell(ic_0,ic_1,2,1)*vectorf(if_0+1,if_1-1)
     &     +prolongcell(ic_0,ic_1,1,2)*vectorf(if_0-1,if_1+1)
     &     +prolongcell(ic_0,ic_1,2,2)*vectorf(if_0+1,if_1+1)
c         if (debug_on.eq.1) then
c           print *, "vectorf[",if_0,",",if_1,"] = ",vectorf(if_0,if_1)
c           call flush(6)
c         endif
        enddo
      enddo

c     if (debug_on.eq.1) then
c       print *, "leaving prolong_am"
c       print *, "vectorf = "
c       call rbugcorn(ff_0,ff_1,lf_0,lf_1,
c    &    ff_0,ff_1,lf_0,lf_1,vectorf)
c       call flush(6)
c     endif
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

      subroutine restrict_am(fc_0,fc_1,lc_0,lc_1,
     &  ff_0,ff_1,lf_0,lf_1,
     &  prolongcell,prolongside0,prolongside1,vectorf, vectorc)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c***********************************************************************
c  Copyright 2006 John A. Trangenstein
c
c  This software is made available for research and instructional use 
c  only. 
c  You may copy and use this software without charge for these 
c  non-commercial purposes, provided that the copyright notice and 
c  associated text is reproduced on all copies.  
c  For all other uses (including distribution of modified versions), 
c  please contact the author at
c    John A. Trangenstein
c    Department of Mathematics
c    Duke University
c    Durham, NC 27708-0320
c    USA
c  or
c    johnt@math.duke.edu
c  
c  This software is made available "as is" without any assurance that it
c  is completely correct, or that it will work for your purposes.  
c  Use the software at your own risk.
c***********************************************************************
      double precision roundoff,small,huge,undefind
      common/machine/roundoff,small,huge,undefind

      double precision e,pi,rt2
      double precision half,one,two,zero
      integer debug_on
      common/const/e,pi,rt2, half,one,two,zero,debug_on
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer fc_0,fc_1,lc_0,lc_1,
     &  ff_0,ff_1,lf_0,lf_1
      double precision 
     &  prolongcell(fc_0:lc_0,fc_1:lc_1,2,2),
     &  prolongside0(fc_0:lc_0+1,fc_1:lc_1,2),
     &  prolongside1(fc_1:lc_1+1,fc_0:lc_0,2),
     &  vectorf(ff_0:lf_0+1,ff_1:lf_1+1)
      double precision 
     &  vectorc(fc_0:lc_0+1,fc_1:lc_1+1)
      integer ic_0,ic_1,if_0,if_1
c     integer i,j
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     if (debug_on.eq.1) then
c       print *, "in restrict_am"
c       print *, "fc = ",fc_0,fc_1
c       print *, "lc = ",lc_0,lc_1
c       print *, "ff = ",ff_0,ff_1
c       print *, "lf = ",lf_0,lf_1
c       print *, "vectorf = "
c       call rbugcorn(ff_0,ff_1,lf_0,lf_1,
c    &    ff_0,ff_1,lf_0,lf_1,vectorf)
c       do j=1,2
c         do i=1,2
c           print *, "prolongcell(",i,",",j,") = "
c           call rbugcell(fc_0,fc_1,lc_0,lc_1,
c    &        fc_0,fc_1,lc_0,lc_1,prolongcell(fc_0,fc_1,i,j))
c         enddo
c       enddo
c       do i=1,2
c         print *, "prolongside0(*,",i,") = "
c         call rbugside0(fc_0,fc_1,lc_0,lc_1,
c    &      fc_0,fc_1,lc_0,lc_1,prolongside0(fc_0,fc_1,i))
c       enddo
c       do i=1,2
c         print *, "prolongside1(*,",i,") = "
c         call rbugside1(fc_0,fc_1,lc_0,lc_1,
c    &      fc_0,fc_1,lc_0,lc_1,prolongside1(fc_1,fc_0,i))
c       enddo
c       call flush(6)
c     endif

c     assumes Dirichlet boundary values
      do ic_1=fc_1+1,lc_1   ! corner index
        if_1=2*ic_1         ! corner index
c       if (debug_on.eq.1) then
c         print *, "ic_1,if_1 = ",ic_1,if_1
c       endif
        do ic_0=fc_0+1,lc_0 ! corner index
          if_0=2*ic_0       ! corner index
c         if (debug_on.eq.1) then
c           print *, "ic_0,if_0 = ",ic_0,if_0
c         endif
          vectorc(ic_0,ic_1)=vectorf(if_0,if_1)
     &      +prolongside0(ic_0,ic_1-1,2)*vectorf(if_0  ,if_1-1)
     &      +prolongside0(ic_0,ic_1  ,1)*vectorf(if_0  ,if_1+1)
     &      +prolongside1(ic_1,ic_0-1,2)*vectorf(if_0-1,if_1  )
     &      +prolongside1(ic_1,ic_0  ,1)*vectorf(if_0+1,if_1  )
     &      +prolongcell(ic_0-1,ic_1-1,2,2)*vectorf(if_0-1,if_1-1)
     &      +prolongcell(ic_0  ,ic_1-1,1,2)*vectorf(if_0+1,if_1-1)
     &      +prolongcell(ic_0-1,ic_1  ,2,1)*vectorf(if_0-1,if_1+1)
     &      +prolongcell(ic_0  ,ic_1  ,1,1)*vectorf(if_0+1,if_1+1)
        enddo
      enddo

c     if (debug_on.eq.1) then
c       print *, "leaving restrict_am"
c       print *, "vectorc = "
c       call rbugcorn(fc_0,fc_1,lc_0,lc_1,
c    &    fc_0,fc_1,lc_0,lc_1,vectorc)
c       call flush(6)
c     endif
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end
