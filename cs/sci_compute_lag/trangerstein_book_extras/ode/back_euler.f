      program main
      integer i,nsteps
      real back_euler,dt,rate,t,tmax,x

      back_euler(x,dt)=x/(1.-dt*rate)

      print *, "enter rate"
      read *, rate
      print *, "enter initial x"
      read *, x
      print *, "enter tmax"
      read *, tmax
      print *, "enter number steps"
      read *, nsteps
      print *, "rate = ",rate," x = ",x
      print *, " tmax = ",tmax," number steps = ",nsteps

      t=0.
      dt=tmax/real(nsteps)
      write(10,*) t,x
      do i=1,nsteps
        x=back_euler(x,dt)
        t=t+dt
        write(10,*) t,x
      enddo

      stop
      end
