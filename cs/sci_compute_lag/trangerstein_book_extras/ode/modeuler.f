      program main
      integer i,nsteps
      real dt,f,rate,t,tmax,u

      f(u,t)=rate*u
      euler(u,t,dt)=u+dt*f(u,t)
      modeuler(u,t,dt)=u+dt*f(euler(u,t,0.5*dt),t+0.5*dt)

      print *, "enter rate"
      read *, rate
      print *, "enter initial u"
      read *, u
      print *, "enter tmax"
      read *, tmax
      print *, "enter number steps"
      read *, nsteps
      print *, "rate = ",rate," u = ",u
      print *, " tmax = ",tmax," number steps = ",nsteps

      t=0.
      dt=tmax/real(nsteps)
      write(12,*) t,u
      do i=1,nsteps
        u=modeuler(u,t,dt)
        t=t+dt
        write(12,*) t,u
      enddo

      stop
      end
