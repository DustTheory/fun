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
c "$Header: /home/faculty/johnt/cvs/flowvar/rcmp.m4,v 1.8 2003/06/15 14:15:28 johnt Exp $"
c "$Header: /home/faculty/johnt/cvs/memdebug/m4amr.sed,v 1.11 2003/06/15 13:50:01 johnt Exp $"
      subroutine rcmpcell(fi_0,fi_1,la_0,la_1,
     &  ibeg_0,ibeg_1,iend_0,iend_1,
     &  array,
     &  arraymax,arraymin)
      integer fi_0,fi_1,la_0,la_1,
     &  ibeg_0,ibeg_1,iend_0,iend_1
      double precision 
     &  array(fi_0:la_0,fi_1:la_1)
      double precision arraymax,arraymin
      double precision roundoff,small,huge,undefind
      common/machine/roundoff,small,huge,undefind

      double precision e,pi,rt2
      double precision half,one,two,zero
      integer debug_on
      common/const/e,pi,rt2, half,one,two,zero,debug_on
      integer ic0,ic1
      integer istart_0,istart_1,ifinish_0,ifinish_1

      arraymax=-huge
      arraymin=huge
      istart_0=max(fi_0,ibeg_0)
      ifinish_0=min(la_0,iend_0)
      istart_1=max(fi_1,ibeg_1)
      ifinish_1=min(la_1,iend_1)
      do ic1=istart_1,ifinish_1
        do ic0=istart_0,ifinish_0
          arraymax=max(arraymax,array(ic0,ic1))
          arraymin=min(arraymin,array(ic0,ic1))
        enddo
      enddo
      return
      end
