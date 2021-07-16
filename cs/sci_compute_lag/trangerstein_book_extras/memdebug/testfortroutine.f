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
c "$Header: /home/faculty/johnt/cvs/deal_new/memdebug/testfortroutine.f,v 1.1 2009/08/20 17:33:34 johnt Exp $"
        subroutine fort(n, a)
        integer n
        double precision a(n)
        external mem_check
        
c       uncommenting the next line will cause a write out of bounds:
c         Fortran starts array addresses at 1, while C/C++ start at 0
c       a(0)=0.d0

        a(n)=0.d0

c       if n is out of bounds for the array declaration in C++, then
c         the next line will cause an error
        call mem_check()

c       it is expensive to call mem_check unnecessarily
c       only use mem_check to find bugs

        return
        end
