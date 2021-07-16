      program main
      integer i,n
      double precision den,fact,mult,pow,quad,rho

      print *, "enter number of iterations"
      read *, n
      print *, "enter reduction factor"
      read *, rho
      print *, "enter multiplier"
      read *, mult

      fact=mult
      pow=mult
      quad=rho
      do i=1,n
        den=1.d0/dble(i)
        fact=fact*den
        pow=pow*rho
        quad=quad*quad
        write(10,10) n,1.d0+mult*den,1.d0+pow,1.d0+fact*pow/mult,
     &    1.d0+mult*quad
   10   format(i2,4d24.16)
      enddo

      write(11,*) "convergent sequence"
      do i=1,n
        write(11,11) dlog(mult/dble(i)),dlog(mult/dble(i+1))
   11   format(2d24.16)
      enddo

      write(12,*) "linearly convergent sequence"
      pow=mult
      do i=1,n
        pow=pow*rho
        write(12,11) dlog(pow),dlog(pow*rho)
      enddo

      write(13,*) "super-linearly convergent sequence"
      fact=mult
      pow=mult
      do i=1,n
        den=1.d0/dble(i)
        fact=fact*den
        pow=pow*rho
        write(13,11) dlog(fact*pow/mult),
     &    dlog(fact*pow*rho/(mult*dble(i+1)))
      enddo

      write(14,*) "quadratically convergent sequence"
      quad=rho
      do i=1,n
        quad=quad*quad
        write(14,11) dlog(mult*quad),dlog(mult*quad*quad)
      enddo

      stop
      end
