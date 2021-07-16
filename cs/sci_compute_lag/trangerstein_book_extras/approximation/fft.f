c     this routine was taken from Cooley, Lewis and Welch, 
c     IEEE transactions E-12 no 1 (March 1965)
c     then had the goto's replaced
      subroutine fft(m,c)
      integer i,ii,j,k,l,m,ndiv2,n,p,twol,twolm1
      double complex c(0:0)
      double complex temp

      j=1
      ndiv2=2**(m-1)
      n=2*ndiv2
c     apply permutation Q:
      do i=1,n-1
        if (i.lt.j) then
          temp=c(j)
          c(j)=c(i)
          c(i)=temp
        endif
        k=ndiv2
        do while (k.lt.j)
          j=j-k
          k=k/2
        enddo
        j=j+k
      enddo

      do l=1,m
        twolm1=2**(l-1)
        twol=2*twolm1
        u=cmplx(1.d0,0.d0)
        theta=pi/twolm1
        z=cmplx(cos(theta),sin(theta))
        do j=1,twolm1
c         apply permutation P:
          do i=j,n,twol
            p=i+twolm1
            temp=c(p)*u
            c(p)=c(i)-temp
            c(i)=c(i)+temp
          enddo
          u=u*z
        enddo
      enddo

      return
      end
