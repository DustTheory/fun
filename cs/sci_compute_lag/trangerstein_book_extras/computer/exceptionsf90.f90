      program main
      real(kind=4) infinity,nanq,nans
      real(kind=8) d,dh,dinfinity,dnanq,dnans,dzero,dmone,done,dtwo
      data infinity /z'7F800000'/
      data nanq /z'7Fc00000'/
      data nans /z'7F800001'/
      data dinfinity /z'7FF0000000000000'/
      data dnanq /z'7FF8000000000000'/
      data dnans /z'7FFc000000000000'/

      dzero=0.d0
      done=1.d0
      dmone=-1.d0
      dtwo=2.d0

      d=done
      dh=dinfinity
      print '("dinfinity = ",g24.16,"   ",z24)',dh,dh

      d=dnans
      print '("signaling nan = ",g24.16,"   ",z24)',d,d

      d=dnanq
      print '("quiet nan = ",g24.16,"   ",z24)',d,d

      d=done/dinfinity
      print '("1.d0/dinfinity = ",g24.16,"   ",z24)',d,d

      d=huge(d)*dtwo
      print '("huge*2.d0 = ",g24.16,"   ",z24)',d,d

      d=dinfinity*dinfinity
      print '("dinfinity*dinfinity = ",g24.16,"   ",z24)',d,d

      d=done/dzero
      print '("1.d0/0.d0 = "g24.16,"   ",z24)',d,d

      d=dzero/dzero
      print '("0.d0/0.d0 = ",g24.16,"   ",z24)',d,d

      d=dinfinity+dinfinity
      print '("dinfinity+dinfinity = ",g24.16,"   ",z24)',d,d

      d=dinfinity-dinfinity
      print '("dinfinity-dinfinity = ",g24.16,"   ",z24)',d,d

      d=dinfinity/dinfinity
      print '("dinfinity/dinfinity = ",g24.16,"   ",z24)',d,d

      d=dinfinity/dzero
      print '("dinfinity/0.d0 = ",g24.16,"   ",z24)',d,d

      d=log(dmone)
      print '("log(-1.d0) = ",g24.16,"   ",z24)',d,d

      d=log(dzero)
      print '("log(0.d0) = ",g24.16,"   ",z24)',d,d

      d=sqrt(dmone)
      print '("sqrt(-1.) = ",g24.16,"   ",z24)',d,d

      d=acos(dtwo)
      print '("arc cos(2.) = ",g24.16,"   ",z24)',d,d

      stop
      end
