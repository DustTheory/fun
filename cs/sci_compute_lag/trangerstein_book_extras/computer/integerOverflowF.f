      program main
      integer*1 i
      integer*2 ii
      integer*4 iii
      integer*8 iiii
      integer*16 iiiii
      integer j

      i=-1
      write(6,10) i,i
   10 format(i3,z3)
      write(6,10) huge(i),huge(i)
      print *, " "

      ii=-1
      write(6,20) ii,ii
   20 format(i6,z6)
      write(6,20) huge(ii),huge(ii)
      print *, " "

      iii=-1
      write(6,40) iii,iii
   40 format(i12,z12)
      write(6,40) huge(iii),huge(iii)
      print *, " "

      iiii=-1
      write(6,80) iiii,iiii
   80 format(i20,z20)
      write(6,80) huge(iiii),huge(iiii)
      print *, " "

      iiiii=-1
      write(6,160) iiiii,iiiii
  160 format(i40,z40)
      write(6,160) huge(iiiii),huge(iiiii)
      print *, " "

      print *, "compute j_n = 2*j_{n-1} + 1, j_0 = 1 while j_n > 0"
      j=1
      do while (j.gt.0)
c       print *,j
        write(6,40) j,j
        j=2*j+1
      enddo
      write(6,40) j,j

      print *
      print *, "compute j_n = 2*j_{n-1} - 1, j_0 = -1 while j_n < 0"
      j=-1
      do while (j.lt.0)
        write(6,40) j,j
        j=2*j-1
      enddo
      write(6,40) j,j

      stop
      end
