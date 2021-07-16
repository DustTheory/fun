      integer function idmax(n,a,inca)
      integer n,inca
      real*8 a(inca,1)
      real*8 amx
      integer i

      idmax=1
      amx=a(1,1)
      do i = 2,n
        if (a(1,i).gt.amx) then
          idmax=i
          amx=a(1,i)
        endif
      enddo
      return
      end
