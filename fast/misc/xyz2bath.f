c
      parameter(nxmax=10001, nymax=10001)
      real interface(nxmax,nymax)
      character*80 filei
c
      open(35, file='for.header', status='old')
c
c     read in model dimensions
c
      read(35,115) xmin,xmax,ymin,ymax,zmin,zmax,size,nx,ny,nz
115   format(7f10.3,3i10)
c
      write(6,1)
1     format('Enter xyz ascii file name')
      read(5,*) filei
c
      open(11, file=filei, status='old')
      open(12, file='int.out', form='unformatted')
c
      nt=nx*ny
c
      do i=1,nt
         read(11,*) x,y,z
         ix=nint((x-xmin)/size)+1
         jy=nint((y-ymin)/size)+1
         interface(ix,jy)=z*1000.
      enddo
c
      do 31 j=1,ny
31          write(12) (interface(i,j),i=1,nx)
c
      stop
      end 
