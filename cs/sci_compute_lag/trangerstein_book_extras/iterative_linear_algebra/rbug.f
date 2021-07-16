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
c "$Header: /home/faculty/johnt/cvs/flowvar/rbug.m4,v 1.7 2004/06/23 21:50:58 johnt Exp $"
c "$Header: /home/faculty/johnt/cvs/memdebug/m4amr.sed,v 1.11 2003/06/15 13:50:01 johnt Exp $"
      subroutine rbugcell(fi_0,fi_1,la_0,la_1,
     &  ibeg_0,ibeg_1,iend_0,iend_1,
     &  array)
      integer fi_0,fi_1,la_0,la_1,
     &  ibeg_0,ibeg_1,iend_0,iend_1
      double precision 
     &  array(fi_0:la_0,fi_1:la_1)
      integer ic0,ic1
      do ic1=ibeg_1,iend_1
        do ic0=ibeg_0,iend_0
          call rwrite_array_2d("array",ic0,ic1,array(ic0,ic1))
        enddo
      enddo
      call flush(6)
      return
      end
