c***********************************************************************
c  Copyright 2006 John A. Trangenstein
c
c  This software is made available for research and instructional use 
c  only. 
c  You may copy and use this software without charge for these 
c  non-commercial purposes, provided that the copyright notice and 
c  associated text is reproduced on all copies.  
c  For all other uses (including distribution of modified versions), 
c  please contact the author at
c    John A. Trangenstein
c    Department of Mathematics
c    Duke University
c    Durham, NC 27708-0320
c    USA
c  or
c    johnt@math.duke.edu
c  
c  This software is made available "as is" without any assurance that it
c  is completely correct, or that it will work for your purposes.  
c  Use the software at your own risk.
c***********************************************************************
c "$Header: /home/faculty/johnt/cvs/flowvar/fortio.m4,v 1.2 2004/06/23 21:50:58 johnt Exp $"
c "$Header: /home/faculty/johnt/cvs/memdebug/m4amr.sed,v 1.11 2003/06/15 13:50:01 johnt Exp $"







      subroutine my_abort()
      double precision bomb
      bomb=0.d0
      bomb=1.d0/bomb
      return
      end
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine write_label(name)
      character*(*) name
      print *, "      ",name
      return
      end
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine nwrite_dimensions(name,a_0,a_1)
      character*(*) name
      integer a_0,a_1
      print *, "      ",name," = ",a_0," ",a_1
      return
      end
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine write_dimensions(name,array)
      character*(*) name
      integer array(0:2-1)
      integer i
      print *, "      ",name," = ",(array(i)," ",i=0,2-1)
      return
      end
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine bwrite_val(name,val)
      character*(*) name
      logical*1 val
      print *,"      ",name," = ",val 
      return
      end
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine bwrite_array_1d(name,ic,array_val)
      character*(*) name
      integer ic
      logical*1 array_val
      print *, "      ",name,"[",ic,"] = ",array_val
      return
      end
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine bwrite_array_2d(name,ic0,ic1,array_val)
      character*(*) name
      integer ic0,ic1
      logical*1 array_val
      print *, "      ",name,"[",ic0,",",ic1,"] = ",array_val
      return
      end
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine bwrite_array_3d(name,ic0,ic1,ic2,array_val)
      character*(*) name
      integer ic0,ic1,ic2
      logical*1 array_val
      print *, "      ",name,"[",ic0,",",ic1,",",ic2,"] = ",array_val
      return
      end
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine bwrite_outs_array_1d(name,array_vall,array_valr)
      character*(*) name
      logical*1 array_vall,array_valr
      print *, "      ",name,"[,1-2] = ",array_vall," ",array_valr
      return
      end
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine bwrite_outs_array_2d(name,ic,array_vall,array_valr)
      character*(*) name
      integer ic
      logical*1 array_vall,array_valr
      print *, "      ",name,"[",ic,",1-2] = ",array_vall," ",array_valr
      return
      end
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine bwrite_outs_array_3d(name,ic1,ic2,array_vall,
     &  array_valr)
      character*(*) name
      integer ic1,ic2
      logical*1 array_vall,array_valr
      print *, "      ",name,"[",ic1,",",ic2,",1-2] = ",array_vall," ",
     &  array_valr
      return
      end
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine iwrite_val(name,val)
      character*(*) name
      integer val
      print *, "      ",name," = ",val
      return
      end
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine iwrite_array_1d(name,ic,array_val)
      character*(*) name
      integer ic,array_val
      print *, "      ",name,"[",ic,"] = ",array_val
      return
      end
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine iwrite_array_2d(name,ic0,ic1,array_val)
      character*(*) name
      integer ic0,ic1,array_val
      print *, "      ",name,"[",ic0,",",ic1,"] = ",array_val
      return
      end
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine iwrite_array_3d(name,ic0,ic1,ic2,array_val)
      character*(*) name
      integer ic0,ic1,ic2,array_val
      print *, "      ",name,"[",ic0,",",ic1,",",ic2,"] = ",array_val
      return
      end
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine iwrite_outs_array_1d(name,array_vall,array_valr)
      character*(*) name
      integer array_vall,array_valr
      print *, "      ",name,"[,1-2] = ",array_vall," ",array_valr
      return
      end
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine iwrite_outs_array_2d(name,ic,array_vall,array_valr)
      character*(*) name
      integer ic,array_vall,array_valr
      print *, "      ",name,"[",ic,",1-2] = ",array_vall," ",array_valr
      return
      end
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine iwrite_outs_array_3d(name,ic1,ic2,
     &  array_vall,array_valr)
      character*(*) name
      integer ic1,ic2,array_vall,array_valr
      print *, "      ",name,"[",ic1,",",ic2,",1-2] = ",array_vall," ",
     &  array_valr
      return
      end

c http://www.math.hawaii.edu/lab/197/fortran/fort3.htm
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine rwrite_val(name,val)
      character*(*) name
      double precision val
c     print *, "      ",name," = ",val
      write(6,10) name,val
   10 format("      ",a," = ",e24.16)
      return
      end
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine rwrite_array_1d(name,ic,array_val)
      character*(*) name
      integer ic
      double precision array_val
c     print *, "      ",name,"[",ic,"] = ",array_val
      write(6,10) name,ic,array_val
   10 format("      ",a,"[",i4,"] = ",e24.16)
      return
      end
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine rwrite_array_2d(name,ic0,ic1,array_val)
      character*(*) name
      integer ic0,ic1
      double precision array_val
c     print *, "      ",name,"[",ic0,",",ic1,"] = ",array_val
      write(6,10) name,ic0,ic1,array_val
   10 format("      ",a,"[",i4,",",i4,"] = ",e24.16)
      return
      end
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine rwrite_array_3d(name,ic0,ic1,ic2,array_val)
      character*(*) name
      integer ic0,ic1,ic2
      double precision array_val
c     print *, "      ",name,"[",ic0,",",ic1,",",ic2,"] = ",array_val
      write(6,10) name,ic0,ic1,ic2,array_val
   10 format("      ",a,"[",i4,",",i4,",",i4,"] = ",e24.16)
      return
      end
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine rwrite_outs_array_1d(name,array_vall,array_valr)
      character*(*) name
      double precision array_vall,array_valr
c     print *, "      ",name,"[,1-2] = ",array_vall," ",array_valr
      write(6,10) name,array_vall,array_valr
   10 format("      ",a,"[,1-2] = ",e24.16," ",e24.16)
      return
      end
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine rwrite_outs_array_2d(name,ic,array_vall,array_valr)
      character*(*) name
      integer ic
      double precision array_vall,array_valr
c     print *, "      ",name,"[",ic,",1-2] = ",array_vall," ",array_valr
      write(6,10) name,ic,array_vall,array_valr
   10 format("      ",a,"[",i4,",1-2] = ",e24.16," ",e24.16)
      return
      end
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine rwrite_outs_array_3d(name,ic1,ic2,
     &  array_vall,array_valr)
      character*(*) name
      integer ic1,ic2
      double precision array_vall,array_valr
c     print *, "      ",name,"[",ic1,",",ic2,",1-2] = ",array_vall," ",
c    &  array_valr
      write(6,10) name,ic1,ic2,array_vall,array_valr
   10 format("      ",a,"[",i4,",",i4,",1-2] = ",e24.16," ",e24.16)
      return
      end
