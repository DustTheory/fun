      subroutine adams_bashforth_startup()
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include 'multistep.i'
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer i,j
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "entering adams_bashforth_startup"
      do i=0,order-1
        lm_gamma(i)=1.d0
        do j=0,i-1
          lm_gamma(i)=lm_gamma(i)-lm_gamma(j)/dble(i-j+1)
        enddo
      enddo
c     print *, "leaving adams_bashforth_startup"
c     print *, "lm_gamma = ",(lm_gamma(i),i=0,order-1)
c     call flush(6)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

      subroutine adams_bashforth(dt,t,y,f,ftable)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include 'multistep.i'
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      double precision dt,t,y
      double precision f
      external f
      double precision ftable(0:order-1)
      double precision fnew,sumf,temp
      integer i
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "entering adams_bashforth"
c     print *, "order = ",order
c     print *, "step = ",step
c     print *, "mem_check 1"
c     call flush(6)
c     call mem_check()
      fnew=f(t,y)
      temp=ftable(order-1)
      ftable(order-1)=fnew
      sumf=fnew
      do i=1,min(step,order-1)
        fnew=ftable(order-i)-temp
        temp=ftable(order-1-i)
        ftable(order-1-i)=fnew
        sumf=sumf+lm_gamma(i)*fnew;
      enddo
      y=y+dt*sumf
      step=step+1
c     print *, "mem_check 2"
c     call flush(6)
c     call mem_check()
c     print *, "leaving adams_bashforth"
c     print *, "ftable = ",(ftable(i),i=order-1,0,-1)
      call flush(6)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

      subroutine adams_moulton_startup()
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include 'multistep.i'
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer i,j
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "entering adams_moulton_startup"
      lm_gamma(0)=1.d0
      do i=1,order-1
        lm_gamma(i)=0.d0
        do j=0,i-1
          lm_gamma(i)=lm_gamma(i)-lm_gamma(j)/dble(i-j+1)
        enddo
      enddo
c     print *, "leaving adams_moulton_startup"
c     print *, "lm_gamma = ",(lm_gamma(i),i=0,order-1)
c     call flush(6)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

      subroutine adams_moulton(dt,t,y,dfdy,f,ftable,ftable_old)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include 'multistep.i'
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      double precision dt,t,y
      double precision dfdy,f
      external dfdy,f
      double precision ftable(0:order-1),ftable_old(0:order-1)
      double precision dfmoultondy,fmoulton,fnew,sumf,sum_gamma,temp,
     &  ynew
      integer i,it
      logical converged
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "entering adams_moulton"
c     print *, "order = ",order
c     print *, "step = ",step
c     print *, "maxit = ",maxit
c     print *, "tolerance = ",tolerance
c     print *, "y = ",y
c     print *, "ftable = ",(ftable(i),i=order-1,0,-1)
c     print *, "mem_check 1"
c     call flush(6)
c     call mem_check()

      sum_gamma=0.d0
      do i=0,min(step,order-1)
        ftable_old(order-1-i)=ftable(order-1-i)
        sum_gamma=sum_gamma+lm_gamma(i)
      enddo
      sum_gamma=sum_gamma*dt
c     print *, "sum_gamma = ",sum_gamma

      converged=.false.
      ynew=y
      it=0
      do while (.not.converged .and. it.le.maxit)
        fnew=f(t+dt,ynew)
c       print *, "f(",ynew,") = ",fnew
        temp=ftable_old(order-1)
        ftable(order-1)=fnew
        sumf=fnew
        do i=1,min(step,order-1)
          fnew=ftable(order-i)-temp
          temp=ftable_old(order-1-i)
          ftable(order-1-i)=fnew
          sumf=sumf+lm_gamma(i)*ftable(order-1-i)
        enddo
c       print *, "ftable = ",(ftable(i),i=order-1,0,-1)
c       print *, "sumf = ",sumf
        fmoulton=ynew-y-dt*sumf
        converged=abs(fmoulton).le.tolerance
     &                           *max(abs(ynew),dt*abs(ftable(order-1)))
        dfmoultondy=1.d0-sum_gamma*dfdy(t+dt,ynew)
        ynew=ynew-fmoulton/dfmoultondy
        it=it+1
c       print *, " "
c       print *, "it = ",it
c       print *, "ynew = ",ynew
c       print *, "fmoulton = ",fmoulton
c       print *, "converged = ",converged
      enddo
      y=ynew
      step=step+1

c     print *, "mem_check 2"
c     call flush(6)
c     call mem_check()
c     print *, "leaving adams_moulton"
c     print *, "y = ",y
c     print *, "ftable = ",(ftable(i),i=order-1,0,-1)
c     call flush(6)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

      subroutine bdf_startup()
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include 'multistep.i'
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer i
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "entering bdf_startup"
      lm_gamma(0)=0.d0
      do i=1,order
        lm_gamma(i)=1./dble(i)
      enddo
c     print *, "leaving bdf_startup"
c     print *, "lm_gamma = ",(lm_gamma(i),i=0,order)
c     call flush(6)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

      subroutine bdf(dt,t,y,dfdy,f,ytable,ytable_old)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include 'multistep.i'
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      double precision dt,t,y
      double precision dfdy,f
      external dfdy,f
      double precision ytable(0:order),ytable_old(0:order)
      double precision dfbdfdy,fbdf,fnew,sumy,sum_gamma,temp,ynew,ytnew
      integer i,it
      logical converged
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "entering bdf"
c     print *, "step = ",step
c     print *, "order = ",order
c     print *, "maxit = ",maxit
c     print *, "tolerance = ",tolerance
c     print *, "ytable = ",(ytable(i),i=order,0,-1)
c     print *, "mem_check 1"
c     call flush(6)
c     call mem_check()

      sum_gamma=0.d0
      do i=0,order
        ytable_old(i)=ytable(i)
        sum_gamma=sum_gamma+lm_gamma(i)
      enddo
      ytable_old(order)=y
c     print *, "sum_gamma = ",sum_gamma

      step=step+1
      converged=.false.
      ynew=y
      it=0
      do while (.not.converged .and. it.le.maxit)
        fnew=f(t+dt,ynew)
c       print *, "f(",ynew,") = ",fnew
        temp=ytable_old(order)
        ytable(order)=ynew
        sumy=0.d0
        do i=1,min(step,order)
          ytnew=ytable(order+1-i)-temp
          temp=ytable_old(order-i)
          ytable(order-i)=ytnew
          sumy=sumy+lm_gamma(i)*ytable(order-i)
        enddo
c       print *, "ytable = ",(ytable(i),i=order,0,-1)
c       print *, "sumy = ",sumy
        fbdf=sumy-dt*fnew
        converged=abs(fbdf).le.tolerance*max(abs(ynew),dt*abs(fnew))
        dfbdfdy=sum_gamma-dt*dfdy(t+dt,ynew)
        ynew=ynew-fbdf/dfbdfdy
        it=it+1
c       print *, " "
c       print *, "it = ",it
c       print *, "ynew = ",ynew
c       print *, "fbdf = ",fbdf
c       print *, "converged = ",converged
      enddo
      y=ynew

c     print *, "mem_check 2"
c     call flush(6)
c     call mem_check()
c     print *, "leaving bdf"
c     print *, "ytable = ",(ytable(i),i=order,0,-1)
c     call flush(6)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

      subroutine explicit_midpoint(dt,t,y,yold,f)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include 'multistep.i'
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      double precision dt,t,y,yold
      double precision f
      external f
      double precision ynew
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     print *, "entering explicit_midpoint"
c     print *, "step = ",step
c     print *, "mem_check 1"
c     call flush(6)
c     call mem_check()
      ynew=yold+2.d0*dt*f(t,y)
      yold=y
      y=ynew
      step=step+1
c     print *, "mem_check 2"
c     call flush(6)
c     call mem_check()
c     print *, "leaving explicit_midpoint"
c     print *, "y,yold = ",y,yold
c     call flush(6)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end
