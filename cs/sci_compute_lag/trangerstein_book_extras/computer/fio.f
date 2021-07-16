      program main
      character*26 capstr
      character*26 str
      logical l
      integer*1 i1
      integer*2 i2
      integer*4 i4
      integer*8 i8
      integer*16 i16
      real*4 r4
      real*8 r8
      integer ii

      character*26 capstr2
      character*26 str2
      character*128 s2
      character*128 s3
      logical ll
      integer*1 ii1
      integer*2 ii2
      integer*4 ii4
      integer*8 ii8
      integer*16 ii16
      real*4 rr4
      real*8 rr8

      l=.false.
      i1=1
      i2=2
      i4=4
      i8=8
      i16=16
      r4=4.
      r8=8.d0
      str='abcdefghijklmnopqrstuvwxyz'
      capstr='00000000000000000000000000'

c     formatted to terminal
      print *, "in fio"
      do ii=1,26
        capstr(ii:ii)=char(ii+64) 
      enddo
      print *, capstr
      print *, str
      print *, "logical = ",l
      print *, "integer*1 = ",i1
      print *, "integer*2 = ",i2
      print *, "integer*4 = ",i4
      print *, "integer*8 = ",i8
      print *, "integer*16 = ",i16
      print *, "real*4 = ",r4
      print *, "real*8 = ",r8

c     formatted to file
      print *, " "
      print *, "formatted write to file"
      open(unit=10,file='fio_formatted_output',status='unknown')
      write(10,*) capstr
      write(10,*) str
      write(10,*) "logical = ",l
      write(10,*) "integer*1 = ",i1
      write(10,*) "integer*2 = ",i2
      write(10,*) "integer*4 = ",i4
      write(10,*) "integer*8 = ",i8
      write(10,*) "integer*16 = ",i16
      write(10,*) "real*4 = ",r4
      write(10,*) "real*8 = ",r8
      close(unit=10,status='keep')

c     formatted from file
      print *, " "
      print *, "formatted read from file"
      open(unit=11,file='fio_formatted_output',status='old');
      read(11,*) capstr2
      print *, capstr2
      read(11,*) str2
      print *, str2
      read(11,*) s2,s3,ll
      print *, "logical = ",ll
      read(11,*) s2,s3,ii1
      print *, "integer*1 = ",ii1
      read(11,*) s2,s3,ii2
      print *, "integer*2 = ",ii2
      read(11,*) s2,s3,ii4
      print *, "integer*4 = ",ii4
      read(11,*) s2,s3,ii8
      print *, "integer*8 = ",ii8
      read(11,*) s2,s3,ii16
      print *, "integer*16 = ",ii16
      read(11,*) s2,s3,rr4
      print *, "real*4 = ",rr4
      read(11,*) s2,s3,rr8
      print *, "real*8 = ",rr8
      close(unit=11,status='keep')

c     unformatted to file
      print *, " "
      print *, "unformatted write to file"
      open(unit=12,file='fio_unformatted_output',form='UNFORMATTED',
     &  status='unknown')
      write(unit=12) capstr
      write(unit=12) str
      write(unit=12) l
      write(unit=12) i1
      write(unit=12) i2
      write(unit=12) i4
      write(unit=12) i8
      write(unit=12) i16
      write(unit=12) r4
      write(unit=12) r8
      close(unit=12,status='keep')

c     unformatted from file
      print *, " "
      print *, "unformatted read from file"
      open(unit=13,file='fio_unformatted_output',form='UNFORMATTED',
     &  status='old');
      read(13) capstr2
      print *, capstr2
      read(13) str2
      print *, str2
      read(13) ll
      print *, "logical = ",ll
      read(13) i1
      print *, "integer*1 = ",ii1
      read(13) i2
      print *, "integer*2 = ",ii2
      read(13) i4
      print *, "integer*4 = ",ii4
      read(13) i8
      print *, "integer*8 = ",ii8
      read(13) i16
      print *, "integer*16 = ",ii16
      read(13) r4
      print *, "real*4 = ",rr4
      read(13) r8
      print *, "real*8 = ",rr8
      close(unit=13,status='keep')

      stop
      end
