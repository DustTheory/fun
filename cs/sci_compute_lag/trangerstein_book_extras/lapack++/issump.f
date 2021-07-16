      integer function issump(n,a,inca)
      integer n,inca
      real a(inca,1)
      integer i
     
      issump=0
      do i = 1,n 
        if (a(1,i).ge.0.d0) issump=issump+1
      enddo   
      return
      end

