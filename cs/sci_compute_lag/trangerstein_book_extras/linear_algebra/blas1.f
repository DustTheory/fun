c find index of first entry of a with max absolute value
      integer function isamax_simple(m,a)
      integer i,m
      real a(m),amx
      isamax_simple=1
      amx=abs(a(1))
      do i=2,m
        if (abs(a(i)).gt.amx) then
          isamax_simple=i
          amx=abs(a(i))
        endif
      enddo
      return
      end

c sum of absolute values of entries of a
      real function sasum_simple(m,a)
      integer i,m
      real a(m)
      do i=1,m
        sasum_simple=sasum_simple+abs(a(i))
      enddo
      return
      end

c y = y + x * alpha
      subroutine saxpy_simple(m,alpha,x,y)
      integer i,m
      real alpha,x(m),y(m)
      do i=1,m
        y(i)=y(i)+x(i)*alpha
      enddo
      return
      end

c y = x
      subroutine scopy_simple(m,x,y)
      integer i,m
      real x(m),y(m)
      do i=1,m
        y(i)=x(i)
      enddo
      return
      end

c return inner product of a and b
      real function sdot_simple(m,a,b)
      integer i,m
      real a(m),b(m)
      sdot_simple=0.
      do i=1,m
        sdot_simple=sdot_simple+a(i)*b(i)
      enddo
      return
      end

c sets all entries of a to alpha
      subroutine sfill(m,alpha,a)
      integer i,m
      real alpha,a(m)
      do i=1,m
        a(i)=alpha
      enddo
      return
      end

c square root of sum of squares of entries of a
      real function snrm2_simple(m,a)
      integer i,m
      real a(m)
      snrm2_simple=0.
      do i=1,m
        snrm2_simple=snrm2_simple+a(i)**2
      enddo
      snrm2_simple=sqrt(snrm2_simple)
      return
      end

c x = x * alpha
      subroutine sscal_simple(m,alpha,x)
      real alpha,x(m)
      integer i,m
      do i=1,m
        x(i)=x(i)*alpha
      enddo
      return
      end

c x <--> y
      subroutine sswap_simple(m,x,y)
      real temp,x(m),y(m)
      integer i,m
      do i=1,m
        temp=x(i)
        x(i)=y(i)
        y(i)=temp
      enddo
      return
      end
