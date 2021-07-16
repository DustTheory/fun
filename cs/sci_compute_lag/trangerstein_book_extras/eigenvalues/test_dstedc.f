      program main
      integer i,info,iwork(28) ! 3+5*n
      double precision d(5),e(4),work(46),Z(5,5) ! 1+n*(4+n)

      do i=1,5
        d(i)=dble(7-i)
      enddo
      do i=1,4
        e(i)=-1.d0
      enddo
      call dstedc('I',5,d,e,Z,5,work,90,iwork,57,info)
      stop
      end
