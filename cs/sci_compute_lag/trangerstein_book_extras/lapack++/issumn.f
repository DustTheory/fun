      integer function issumn(n,a,inca)
      integer n,inca
      real a(inca,1)
      integer i
     
      issumn=0
      do i = 1,n 
        if (a(1,i).le.0.d0) issumn=issumn+1
      enddo   
      return  
      end     

