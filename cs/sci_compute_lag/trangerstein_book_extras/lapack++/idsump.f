      integer function idsump(n,a,inca)
      integer n,inca
      double precision a(inca,1)
      integer i
     
      idsump=0
      do i = 1,n 
        if (a(1,i).ge.0.d0) idsump=idsump+1
      enddo   
      return
      end

