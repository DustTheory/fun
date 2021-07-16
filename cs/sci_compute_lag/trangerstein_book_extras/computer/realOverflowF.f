      program main
      real*4 r4
      real*8 r8
      real a,aold,anew

      r4=-1.
      write(6,40) r4,r4
   40 format(g32.24,2x,z8)
      print *, " "

      r8=-1.
      write(6,80) r8,r8
   80 format(g56.48,2x,z16)
      print *, " "

      print *,"compute a_n = 0.5*a_{n-1} + 1., a_1 = 1",
     &  " until rounding error"
      aold=1.
      anew=1.+0.5*aold
      do while (anew .lt. 2.)
        write(6,40) anew,anew
        anew=1.+0.5*anew
      enddo
      write(6,40) anew,anew

      print *," "
      print *,"compute a_n = 2.*a_{n-1} + 1., a_1 = 1 until overflow"
      a=1.
      do while (.not. (isnan(a) .or. abs(a) > huge(a)))
        write(6,40) a,a
        a=2.*a+1.
      enddo
      write(6,40) a,a

      print *," "
      print *,"compute a_n = 0.5*a_{n-1}, a_1 = 2-epsilon",
     &  " while positive"
      do while (aold.gt.0.)
        write(6,40) aold,aold
        aold=0.5*aold
      enddo
      write(6,40) aold,aold

      stop
      end
