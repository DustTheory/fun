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
c "$Header: /home/faculty/johnt/cvs/deal_new/memdebug/testtypes.m4,v 1.1 2009/08/20 17:33:33 johnt Exp $"

        subroutine test_types()
        double precision a(2)
        real b(2)
        double complex c(2)
        integer d(2)
        logical*1 e(2)

        call check_double_size(a,a(2))
        call check_single_size(b,b(2))
        call check_complex_size(c,c(2))
        call check_integer_size(d,d(2))
        call check_logical_size(e,e(2))

        return
        end
