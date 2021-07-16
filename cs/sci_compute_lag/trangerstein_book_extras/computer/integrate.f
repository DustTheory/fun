      subroutine integrate(nsteps,tmax, soln,time)
c input:
      integer nsteps
      double precision tmax
c input-output:
      double precision soln(0:1),time(0:1)

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

      dt=tmax/float(nsteps)
      t=time(0)
      x=soln(0)
      do i=1,nsteps
        x=method(t,x,dt)
        t=t+dt
        soln(i)=x
        time(i)=t
      enddo

      return
      end
