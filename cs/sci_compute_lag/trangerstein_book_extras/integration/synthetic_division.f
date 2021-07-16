      subroutine synthetic_division(n,coefs,x, p,pderiv)
      integer n
      double precision coefs(0:n),x
      double precision p,pderiv
      integer i

      p=coefs(0)
      pderiv=coefs(0)
      do j=1,n
        p=coefs(j)+x*p
        pderiv=p+x*pderiv
      enddo

      return
      end
