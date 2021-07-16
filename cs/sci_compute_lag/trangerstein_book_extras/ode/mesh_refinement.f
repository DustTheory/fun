      program main
      integer i,j,nsteps,ratio,nlevels
      real dt,error,euler,exact,f,rate,t,tmax,u,uinit

      f(u,t)=rate*u
      euler(u,t,dt)=u+dt*f(u,t)

c try 
c  rate=1.
c  initial u=1.
c  tmax=2.
c  initial number steps = 10
c  refinement ratio = 10
c  number of refinement levels = 7

      print *, "enter rate"
      read *, rate
      print *, "enter initial u"
      read *, uinit
      print *, "enter tmax"
      read *, tmax
      print *, "enter initial number steps"
      read *, nsteps
      print *, "enter refinement ratio"
      read *, ratio
      print *, "enter number of refinement levels"
      read *, nlevels
      print *, "rate = ",rate," uinit = ",uinit
      print *, " tmax = ",tmax," number steps = ",nsteps
      print *, " ratio = ",ratio," nlevels = ",nlevels

      exact=uinit*exp(rate*tmax)
      do j=1,nlevels
        t=0.
        dt=tmax/real(nsteps)
        u=uinit
        do i=1,nsteps
          u=euler(u,t,dt)
          t=t+dt
        enddo
        error=abs(exact-u)
        write(11,*) -log10(dt),log10(error)
        nsteps=nsteps*ratio
      enddo

      stop
      end
