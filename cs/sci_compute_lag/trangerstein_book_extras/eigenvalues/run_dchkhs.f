      program main
      integer i,j,nsizes,nn(1),ntypes,iseed(4),nounit,lda,ldu,
     &  nwork,iwork(100),info
      logical dotype(16),select(100)
      double precision thresh,a(100,100),h(100,100),
     &  t1(100,100),t2(100,100),u(100,100),z(100,100),uz(100,100),
     &  wr1(100),wi1(100),wr3(100),wi3(100),evectl(100,100),
     &  evectr(100,100),evecty(100,100),evectx(100,100),
     &  uu(100,100),tau(100),work(40002),result(14)

      nsizes=1
      nn(1)=100
      ntypes=16
      do i=1,15
        dotype(i)=.false.
      enddo
      dotype(16)=.true.
      do i=1,4
        iseed(i)=i
      enddo
      thresh=16.d0*epsilon(1.d0)
      nounit=10
      lda=100
      ldu=100
      nwork=40002
      call dchkhs(nsizes,nn,ntypes,dotype,iseed,thresh,
     &  nounit,a,lda,h,t1,t2,u,ldu,z,uz,wr1,
     &  wi1,wr3,wi3,evectl,evectr,evecty,evectx,
     &  uu,tau,work,nwork,iwork,select,result,
     &  info)
      do j=1,100
        do i=1,100
          print *, a(i,j)
        enddo
      enddo
      stop
      end
