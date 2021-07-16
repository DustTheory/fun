      function exp1(x)
      double precision exp1,x

c     if (abs(x).le.0.000244140625d0) then
c       exp1=(6.d0+x)/(6.d0-2.d0*x)

c     if (abs(x).le.0.0078125d0) then
c       exp1=(60.d0+x*(6.d0+x))/(60.d0-x*(24.d0-3.d0*x))

c     if (abs(x).le.0.015625d0) then
c       exp1=(840.d0+x*(60.d0+x*(20.d0+x)))
c    &      /(840.d0-x*(360.d0-x*(60.d0-x*4.d0)))

      if (abs(x).le.0.25d0) then
        exp1=(15120.d0+x*( 840.d0+x*( 420.d0+x*( 20.d0+x     ))))
     &      /(15120.d0-x*(6720.d0-x*(1260.d0-x*(120.d0-x*5.d0))))
      else
        exp1=(exp(x)-1.d0)/x
      endif

      return
      end

      function exp1p(x)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      double precision exp1p,x
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     if (abs(x).le.0.000244140625d0) then
c       exp1p=(3.d0+2.d0*x)/(6.d0+3.d0*x)

c     if (abs(x).le.0.00048828125d0) then
c       exp1p=(24.d0+x*(16.d0+x*6.d0))/(48.d0+x*(24.d0+x*8.d0))

c     if (abs(x).le.0.001953125d0) then
c       exp1p=(60.d0+x*(40.d0+x*(15.d0+x*4.d0)))
c    &      /(120.d0+x*(60.d0+x*(20.d0+x*5.d0)))

c     if (abs(x).le.0.0156125d0) then
c       exp1p=(360.d0+x*(240.d0+x*(90.d0+x*(24.d0+x*5.d0))))
c    &      /(720.d0+x*(360.d0+x*(120.d0+x*(30.d0+x*6.d0))))

      if (abs(x).le.0.03125d0) then
        exp1p=
     &    (2520.d0+x*(1680.d0+x*(630.d0+x*(168.d0+x*(35.d0+x*6.d0)))))
     &   /(5040.d0+x*(2520.d0+x*(840.d0+x*(210.d0+x*(42.d0+x*7.d0)))))
      else
        exp1p=(exp(x)*(x-1.d0)+1.d0)/(x*(exp(x)-1.d0))
      endif
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

      function time_independent_potassium_ak1(vm)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include 'time_independent_k.i'
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      double precision time_independent_potassium_ak1,vm
      double precision vmminusek1,argak1
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      vmminusek1=vm-time_independent_k_reversal_potential
      argak1=expmultak1*(vmminusek1+expshiftak1)
      if (argak1.ge.0.d0) then
        time_independent_potassium_ak1=
     &      multak1*exp(-argak1)/(1.d0+exp(-argak1))
      else
        time_independent_potassium_ak1=multak1/(1.d0+exp(argak1))
      endif
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

      function time_independent_potassium_ak1_derivative(vm)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include 'time_independent_k.i'
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      double precision time_independent_potassium_ak1_derivative,vm
      double precision time_independent_potassium_ak1
      external time_independent_potassium_ak1
      double precision vmminusek1,argak1
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      vmminusek1=vm-time_independent_k_reversal_potential
      argak1=expmultak1*(vmminusek1+expshiftak1)
      if (argak1.ge.0.d0) then
        time_independent_potassium_ak1_derivative=
     &    -((time_independent_potassium_ak1(vm)*expmultak1)
     &    /(1.d0+exp(-argak1)))
      else
        time_independent_potassium_ak1_derivative=
     &    -((time_independent_potassium_ak1(vm)*expmultak1)
     &    *exp(argak1)/(1.d0+exp(argak1)))
      endif
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

c***********************************************************************

      function time_independent_potassium_bk1(vm)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include 'time_independent_k.i'
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      double precision time_independent_potassium_bk1,vm
      double precision vmminusek1,argbk1n1,argbk1n2,argbk1d
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      vmminusek1=vm-time_independent_k_reversal_potential
      argbk1n1=expmultbk1n1*(vmminusek1+expshiftbk1n1)
      argbk1n2=expmultbk1n2*(vmminusek1+expshiftbk1n2)
      argbk1d=expmultbk1d*(vmminusek1+expshiftbk1d)
      if (argbk1d.ge.0.d0) then
        time_independent_potassium_bk1=
     &      (multbk1n1*exp(argbk1n1-argbk1d)+exp(argbk1n2-argbk1d))
     &      /(1.d0+exp(-argbk1d))
      else
        time_independent_potassium_bk1=
     &      (multbk1n1*exp(argbk1n1)+exp(argbk1n2))
     &      /(1.d0+exp(argbk1d))
      endif
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end

      function time_independent_potassium_bk1_derivative(vm)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include 'time_independent_k.i'
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      double precision time_independent_potassium_bk1_derivative,vm
      double precision time_independent_potassium_bk1
      external time_independent_potassium_bk1
      double precision vmminusek1,argbk1n1,argbk1n2,argbk1d
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      vmminusek1=vm-time_independent_k_reversal_potential
      argbk1n1=expmultbk1n1*(vmminusek1+expshiftbk1n1)
      argbk1n2=expmultbk1n2*(vmminusek1+expshiftbk1n2)
      argbk1d=expmultbk1d*(vmminusek1+expshiftbk1d)
      if (argbk1d.ge.0.d0) then
        time_independent_potassium_bk1_derivative=
     &   ( (multbk1n1*exp(argbk1n1)*expmultbk1n1
     &     +exp(argbk1n2)*expmultbk1n2)
     &     -time_independent_potassium_bk1(vm)*exp(argbk1d)*expmultbk1d)
     &      /(1.d0+exp(argbk1d))
      else
        time_independent_potassium_bk1_derivative=
     &    ( (multbk1n1*exp(argbk1n1-argbk1d)*(expmultbk1n1-expmultbk1d)
     &      +exp(argbk1n2-argbk1d)*(expmultbk1n2-expmultbk1d))
     &    +time_independent_potassium_bk1(vm)*expmultbk1d*exp(-argbk1d))
     &      /(1.d0+exp(-argbk1d))
      endif
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end
