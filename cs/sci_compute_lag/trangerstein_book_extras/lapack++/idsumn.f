      integer function idsumn(n,a,inca)
      integer n,inca
      double precision a(inca,1)
      integer i
     
      idsumn=0
      do i = 1,n 
        if (a(1,i).le.0.d0) idsumn=idsumn+1
      enddo   
      return  
      end     

