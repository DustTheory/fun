      program main
      real(kind=4) r4
      real(kind=8) r8
      real a,aold,anew

      r4=-1.
      print '(g32.24,2x,z8)',r4,r4
      print '(g32.24,2x,z8)',huge(1._4),huge(1._4)
      print *, " "

      r8=-1.
      print '(g56.48,2x,z16)',r8,r8
      print '(g56.48,2x,z16)',huge(1._8),huge(1._8)
      print *, " "

      print *,"compute a_n = 0.5*a_{n-1} + 1., a_1 = 1", &
     &  " until rounding error"
      aold=1.
      anew=1.+0.5*aold
      do while (anew .lt. 2.)
        print '(g32.24,2x,z8)',anew,anew
        anew=1.+0.5*anew
      enddo
      print '(g32.24,2x,z8)',anew,anew

      print *," "
      print *,"compute a_n = 2.*a_{n-1} + 1., a_1 = 1 until overflow"
      a=1.
      do while (.not. (isnan(a) .or. abs(a) > huge(a)))
        print '(g32.24,2x,z8)',a,a
        a=2.*a+1.
      enddo
      print '(g32.24,2x,z8)',a,a

      print *," "
      print *,"compute a_n = 0.5*a_{n-1}, a_1 = 2-epsilon", &
     &  " while positive"
      do while (aold.gt.0.)
        print '(g32.24,2x,z8)',aold,aold
        aold=0.5*aold
      enddo
      print '(g32.24,2x,z8)',aold,aold

      stop
      end
