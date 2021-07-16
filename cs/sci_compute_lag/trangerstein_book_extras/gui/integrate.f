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
c "$Header: /home/faculty/johnt/cvs/deal_new/gui/integrate.f,v 1.1 2009/08/20 17:32:37 johnt Exp $"
      subroutine integrate(nsteps,tmax, soln,time)
c input:
      integer nsteps
      double precision tmax
c input-output:
      double precision soln(0:nsteps),time(0:nsteps)
c common block
      double precision rate
      common/odepar/rate
c local
      integer i
      double precision dt,method,f,t,x
c arithmetic statement functions
      f(t,x)=rate*x
      method(t,x,dt)=x+dt*f(t,x)

c     print *, "nsteps,tmax = ",nsteps,tmax
c     print *, "soln,time = ",soln(0),time(0)
c     print *, "rate = ",rate
c     call flush(6)

      dt=tmax/dfloat(nsteps)
      t=time(0)
      x=soln(0)
      do i=1,nsteps
        x=method(t,x,dt)
        t=t+dt
        soln(i)=x
        time(i)=t
c       print *, "i,t,x = ",i,t,x
      enddo

      return
      end
