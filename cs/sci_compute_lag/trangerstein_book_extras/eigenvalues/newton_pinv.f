      program main
      integer i,j,its
      double precision alpha,A(3,2),AXj(3),delta,sigma,X(2,3),XAX(2,3)
      double precision beta,rt2,U(3,3),V(2,2),VhXU(2,3),XUj(2)

      delta=sqrt(epsilon(1.d0))
      sigma=sqrt(2.d0+delta*delta)
      alpha=0.75d0
      beta=1.d0/sqrt(2.d0+delta*delta)
      rt2=sqrt(2.d0)

      A(1,1)=1.d0
      A(2,1)=delta
      A(3,1)=0.d0
      A(1,2)=1.d0
      A(2,2)=0.d0
      A(3,2)=delta

      U(1,1)=beta*sqrt(2.d0)
      U(2,1)=delta*beta/rt2
      U(3,1)=U(2,1)
      U(1,2)=0.d0
      U(2,2)=-1.d0/rt2
      U(3,2)=1.d0/rt2
      U(1,3)=-delta*beta
      U(2,3)=beta
      U(3,3)=beta
      V(1,1)=1.d0/rt2
      V(2,1)=V(1,1)
      V(1,2)=-V(1,1)
      V(2,2)=V(1,1)

      do j=1,2
        do i=1,3
          X(j,i)=A(i,j)*alpha
        enddo
      enddo
c     print *, "X = "
c     do i=1,2
c       print *, (X(i,j),j=1,3)
c     enddo

      do its=1,60
        do j=1,3
          call dgemv('N',3,2,1.d0,A,3,X(1,j),1,0.d0,AXj,1)
          call dgemv('N',2,3,1.d0,X,2,AXj,1,0.d0,XAX(1,j),1)
        enddo
        call dscal(6,2.d0,X,1)
        call daxpy(6,-1.d0,XAX,1,X,1)
        print *, "its = ",its," X = "
        do i=1,2
          print *, (X(i,j),j=1,3)
        enddo

        do j=1,3
          call dgemv('N',2,3,1.d0,X,2,U(1,j),1,0.d0,XUj,1)
          call dgemv('T',2,2,1.d0,V,2,XUj,1,0.d0,VhXU(1,j),1)
        enddo
        print *, "V^T X U = "
        do i=1,2
          print *, (VhXU(i,j),j=1,3)
        enddo
      enddo

      stop
      end
