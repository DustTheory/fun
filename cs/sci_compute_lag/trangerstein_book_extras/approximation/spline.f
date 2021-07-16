      subroutine spline_setup
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include 'spline.i'
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      double precision 
     &  time_independent_potassium_ak1,
     &    time_independent_potassium_ak1_derivative,
     &  time_independent_potassium_bk1,
     &    time_independent_potassium_bk1_derivative
      external 
     &  time_independent_potassium_ak1,
     &    time_independent_potassium_ak1_derivative,
     &  time_independent_potassium_bk1,
     &    time_independent_potassium_bk1_derivative
      integer i
      double precision vm
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      dpotential_spline=(potential_max_spline-potential_min_spline)
     &                 /dble(nspline)
      do i=0,nspline
        vm=potential_min_spline+dble(i)*dpotential_spline
        potential_spline(i)=vm

        ak1_spline(i)=time_independent_potassium_ak1(vm)
        ak1_derivative_spline(i)=
     &    time_independent_potassium_ak1_derivative(vm)
     &    *dpotential_spline
        bk1_spline(i)=time_independent_potassium_bk1(vm)
        bk1_derivative_spline(i)=
     &    time_independent_potassium_bk1_derivative(vm)
     &    *dpotential_spline
      enddo
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

c***********************************************************************

      subroutine spline_values(vm, ak1,bk1)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include 'spline.i'
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      double precision vm
      double precision ak1,bk1
      integer i
      double precision onemy,onemy2,s0,s1,v0,v1,y,y2
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      i=int((vm-potential_min_spline)/dpotential_spline)
      y=(vm-potential_spline(i))/dpotential_spline
      y2=y*y
      onemy=1.d0-y
      onemy2=onemy*onemy
      v0=(1.d0+2.d0*y)*onemy2
      s0=y*onemy2
      v1=(3.d0-2.d0*y)*y2
      s1=-(onemy*y2)

      ak1=v0 * ak1_spline(i)   + s0 * ak1_derivative_spline(i)
     &  +v1 * ak1_spline(i+1) + s1 * ak1_derivative_spline(i+1)
      bk1=v0 * bk1_spline(i)   + s0 * bk1_derivative_spline(i)
     &  +v1 * bk1_spline(i+1) + s1 * bk1_derivative_spline(i+1)

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

c***********************************************************************

      subroutine spline_derivative_values(vm, ak1,bk1)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include 'spline.i'
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      double precision vm
      double precision ak1,bk1
      integer i
      double precision onemy,s0,s1,v0,v1,y
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      i=(vm-potential_min_spline)/dpotential_spline
      y=(vm-potential_spline(i))/dpotential_spline
      onemy=1.d0-y
      v0=-(6.d0*y*onemy)
      s0=onemy*(1.d0-3.d0*y)
      v1=-v0
      s1=y*(3.d0*y-2.d0)

      ak1=(v0 * ak1_spline(i)   + s0 * ak1_derivative_spline(i)
     &   +v1 * ak1_spline(i+1) + s1 * ak1_derivative_spline(i+1))
     &  /dpotential_spline
      bk1=(v0 * bk1_spline(i)   + s0 * bk1_derivative_spline(i)
     &   +v1 * bk1_spline(i+1) + s1 * bk1_derivative_spline(i+1))
     &  /dpotential_spline

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end
