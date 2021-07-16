      program main
      integer nbig
      parameter (nbig=5500)
      integer i,j,k,m,n,nonzero
      double precision A(nbig,nbig)

      open (unit=15,file='s1rmq4m1.mtx',status='old')
      read(15,*) m,n,nonzero
      do k=1,nonzero
        read(15,*) i,j,A(i,j)
        A(j,i)=A(i,j)
      enddo
      stop
      end
