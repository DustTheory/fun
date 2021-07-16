      program hw19
      integer bign,m,n,nonzero,i,j,k,LWORK,INFO
      parameter(bign = 3140,LWORK=386220)
      double precision A(bign,bign),S(bign),U(bign,bign)
      double precision VT(bign,bign),WORK(LWORK)
      double precision EYE1(bign,bign),EYE2(bign,bign),NORM1,NORM2,NORM3
      double precision dlange
      double precision Ac(bign,bign)
      real elapsed,reslt,tarray(2),start

      open(unit=191, file = 'psmigr_1.mtx', status='old')
      read(191,*) m,n,nonzero
      print *,'m, n, nonzero = ', m,n,nonzero
c
c
c     FILL MATRIX WITH DATA VALUES
      call dtime(tarray,reslt)
      start=tarray(1)
      do k=1,nonzero
        read(191,*) i, j, A(i,j)
      enddo
      call dtime(tarray,reslt)
      elapsed=tarray(1)-start
      print *, "reading matrix took ",elapsed," seconds"
        
c     COPY MATRIX A
      call dlacpy('A',m,n,A,bign,Ac,bign)

        
c     COMPUTE SVD
c     call dgesvd('A','A',m,n,A,bign,S,U,bign,VT,bign,work,-1,INFO)
c     LWORK=int(work(1))
c     print *, "LWORK = ",LWORK
      call dtime(tarray,reslt)
      start=tarray(1)
      call dgesvd('A','A',m,n,A,bign,S,U,bign,VT,bign,WORK,LWORK,INFO)
      call dtime(tarray,reslt)
      elapsed=tarray(1)-start
      print *, "calling dgesvd took ",elapsed," seconds"
        
c     print *,'INFO = ',INFO
c     print *,'optimal LWORK = ',WORK(1)
c     print *, 'U',U
c     print *, 'S',S
c     print *, 'VT',VT
        
c     FIND F NORM OF U**H*U - I
c     compute U**H*U - I
      call dtime(tarray,reslt)
      start=tarray(1)
      call dgemm('T','N',m,m,m,1.d0,U,bign,U,bign,0.d0,EYE1,bign)
        
      do i = 1,n
        EYE1(i,i) = EYE1(i,i) - 1.d0
      enddo 
        
      NORM1 = dlange('F',m,m,EYE1,bign,WORK)
      print *,'FNORM(U**H*U-I) = ',NORM1
        
c     FIND F NORM OF V**H*V - I
c     compute V**H*V - I
      call dgemm('N','T',n,n,n,1.d0,VT,bign,VT,bign,0.d0,EYE2,bign)
        
      do i = 1,n
        EYE2(i,i) = EYE2(i,i) - 1.d0
      enddo 
                
      NORM2 = dlange('F',n,n,EYE2,bign,WORK)
      print *,'FNORM(V**H*V-I) = ',NORM2
                
c     FIND F NORM OF A - U*S*V**H
      do j=1,n
        call dscal(m,S(j),U(1,j),1)
      enddo
      call dgemm('N','N',m,n,m,-1.d0,U,bign,VT,bign,1.d0,Ac,bign)
        
      NORM3 = dlange('F',m,n,Ac,bign,WORK)
      print *, 'FNORM(A-U*S*V**H) = ',NORM3
      call dtime(tarray,reslt)
      elapsed=tarray(1)-start
      print *, "computing norms took ",elapsed," seconds"
        
      stop
      end

