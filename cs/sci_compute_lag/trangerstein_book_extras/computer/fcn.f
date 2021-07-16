      subroutine fcn1(n,A)
      integer n
      double precision A(n)
      integer i
      double precision temp

      temp=2.2d0**3.3d0
      do i=1,n
        A(i)=temp
      enddo

      return
      end

      subroutine fcn2(n,A)
      integer n
      double precision A(n)
      double precision cubed,x
      integer i
      cubed(x)=x**3.3d0

      do i=1,n
        A(i)=cubed(2.2d0)
      enddo

      return
      end

      subroutine fcn3(A)
      double precision A

      A=2.2d0**3.3d0

      return
      end

      function fcn4(x)
      double precision fcn4,x
      fcn4=x**3.3d0
      return
      end
