      program main
      real(kind=8) a,b,c,m,s,disc,large,small

      print *, "enter a,b,c for quadratic a * x^2 + b * x + c = 0"
      read *, a,b,c
      print *, "a,b,c = ",a,b,c

      m =max( abs( a ), abs( b ), abs( c ) )
      if ( m .le. 0.d0 ) then
        print *, "all coefficients zero, so roots are arbitrary"
      else
        s = 2.d0 ** int( log( m ) / log( 2.d0 ) )
        if ( s < scale ) s = s * 2.d0
        a = a / s
        b = b / s
        c = c / s
        if ( abs( a ) .gt. 0.d0 ) then
          b = b / ( -2.d0 * a )
          c = c / a
          disc = b * b - c
          if ( disc .ge. 0.d0 ) then
            disc = sqrt( disc )
            large = b + sign(disc,b)
            if ( abs(large) .gt. 0.d0 ) then
              small = c / large
            else
              small=0.d0
            endif
            print '("two real roots = ",g24.16," ",g24.16)',large,small
          else
            print '("complex roots: real part = ",g24.16," imaginary part = ",g24.16)',b,sqrt(-disc)
          endif
        else
          if ( abs( b ) .gt. 0.d0 ) then
            print '("one real root = ",g24.16)',c/b
          else
            print *, "no roots: a = 0 = b"
          endif
        endif
      endif

      stop
      end
