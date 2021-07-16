      subroutine gram_schmidt(A,m,n,R)
c     Q is stored in A
c     diagonal of R contains norm squared of columns of Q
      integer j,k,m,n
      real A(m,n),R(n,n)
      do k=1,n
        R(k,k)=sdot(m,A(1,k),1,A(1,k),1)
        if (k.lt.n) then
          call sgemv('T',m,n-k,1./R(k,k),A(1,k+1),m,
     &      A(1,k),1,0.,R(k,k+1),n)
          call sger(m,n-k,-1.,A(1,k),1,R(k,k+1),n,A(1,k+1),m)
        endif
c       do j=k+1,n
c         R(k,j)=sdot(m,A(1,k),1,A(1,j),1)
c         call saxpy(m,-R(k,j),A(1,k),1,A(1,j),1)
c       enddo
      enddo
      return
      end

      subroutine successive_orthogonal_projection(Q,m,n,k,R,B,Y)
c     residuals are stored in B, which corresponds to k right-hand sides
      integer j,k,m,n
      real Q(m,n),R(n,n),B(m,k),Y(n,k)
      do j=1,n
        call sgemv('T',m,k,1./R(j,j),B,m,A(1,j),1,0.,Y(j,1),n)
        call sger(m,k,-1.,A(1,j),1,Y(j,1),n,B,m)
      enddo
c     do j=1,n
c       y(j)=sdot(m,Q(1,j),1,b,1)/R(j,j)
c       call saxpy(m,-y(j),Q(1,j),1,b,1)
c     enddo
      return
      end

c     can solve Rx=y for successive orthogonal projection
c     by using LAPACK BLAS routine strsv as follows:
c       call strsv('U','N','U',n,R,n,y,1)
c     This routine stores x in y

      subroutine classical_gram_schmidt(A,m,n,R)
c     Q is stored in A
      integer j,k,m,n
      real A(m,n),R(n,n)
      do k=1,n
        if (k.gt.1) then
          call sgemv('T',m,k-1,1.,A,m,A(1,k),1,0.,R(1,k),1)
          call sgemv('N',m,k-1,-1.,A,m,R(1,k),1,1.,A(1,k),1)
        endif
c       do j=1,k-1
c         R(j,k)=sdot(m,A(1,j),1,A(1,k),1)
c         call saxpy(m,-R(j,k),A(1,j),1,A(1,k),1)
c       enddo
        R(k,k)=snrm2(m,A(1,k),1)
        call sscal(m,1./R(k,k),A(1,k),1)
      enddo
      return
      end

      subroutine simultaneous_orthogonal_projection(Q,m,n,b,y)
c     residual is stored in b
      integer j,m,n
      real Q(m,n),b(m),y(n)
      call sgemv('T',m,n,1.,Q,m,b,1,0.,y,1)
      call sgemv('N',m,n,-1.,Q,m,y,1,1.,b,1)
c     do j=1,n
c       y(j)=sdot(m,Q(1,j),1,b,1)
c     enddo
c     do j=1,n
c       call saxpy(m,-y(j),Q(1,j),1,b,1)
c     enddo
      return
      end

c     for simultaneous orthogonal projection, can solve Rx=y 
c     by using LAPACK BLAS routine strsv as follows:
c       call strsv('U','N','N',n,R,n,y,1)
c     This routine stores x in y
