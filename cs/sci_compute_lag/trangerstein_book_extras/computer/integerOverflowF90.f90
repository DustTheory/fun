      program main
      integer(kind=1) i
      integer(kind=2) ii
      integer(kind=4) iii
      integer(kind=8) iiii
      integer(kind=16) iiiii
      integer j

      i=-1
      write(6,10) i,i
   10 format(i3,z3)
      write(6,10) huge(1_1),huge(1_1)
      print *, " "

      ii=-1
      write(6,20) ii,ii
   20 format(i6,z6)
      write(6,20) huge(1_2),huge(1_2)
      print *, " "

      iii=-1
      write(6,40) iii,iii
   40 format(i12,z12)
      write(6,40) huge(1_4),huge(1_4)
      print *, " "

      iiii=-1
      write(6,80) iiii,iiii
   80 format(i20,z20)
      write(6,80) huge(1_8),huge(1_8)
      print *, " "

      iiiii=-1
      write(6,160) iiiii,iiiii
  160 format(i40,z40)
      write(6,160) huge(1_16),huge(1_16)
      print *, " "

      print *, "compute j_n = 2*j_{n-1} + 1, j_0 = 1 while j_n > 0"
      j=1
      do while (j.gt.0)
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
