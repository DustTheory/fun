      subroutine f90copy(m,n,A,B)
      integer :: m,n
      real(kind=8), intent(in) :: A(m,n)
      real(kind=8), intent(out) :: B(m,n)
      B=A
      return
      end
