      program main
c declare variables:
      integer i,nsteps
      real dt,f,rate,t,tmax,u

c arithmetic statement functions:
      f(u,t)=rate*u
      euler(t,u,dt)=u+dt*f(u,t)

c read input:
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

c apply forward-Euler method:
      t=0.
      dt=tmax/real(nsteps)
      write(10,*) t,u
      do i=1,nsteps
        u=euler(t,u,dt)
        t=t+dt
        write(10,*) t,u
      enddo

      stop
      end
