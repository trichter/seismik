c
c     version 1.0  Apr 1996
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |                ********  S M O D E L  *****                  |
c     |                                                              |
c     |    Convert a velocity model to a slowness model using the    |
c     |                       inverse grid                           |
c     |                                                              |
c     |                 Written by Colin A. Zelt                     |
c     |                                                              |
c     |                Dept of Geology & Geophysics                  |   
c     |                      Rice University                         |
c     |                  Houston, TX  77005-1892                     |   
c     |                                                              |
c     ----------------------------------------------------------------
c
c        patched for g77 compilation
c
c        Drew Brenders (brenders@geoladm.geol.queensu.ca)
c
c
c     I/O units:
c
c        10 -- input:  inverse 3D data grid
c
c        11 -- output:  forward 3D data grid
c
c        35 -- input:  parameters describing size of forward 3D grid
c
c        36 -- input:  parameters describing size of inverse 3D grid
c
c
c     ----------------------------------------------------------------
c
c
      include 'ray.par'
      common /blk/ xmin,ymin,zmin,xinc,yinc,zinc,nxi,nyi,nzi,veln
c
      real smod(nximax,nyimax,nzimax)
      integer*2 ivel(nxmax,nymax,nzmax)
      integer cellstotal
      character*72 file1,file2
c
      write(6,335)
335   format('SMODEL: convert velocity model to slowness model')
      open(50, file='nowrite', status='old', err=9999, iostat=ioflag)
9999  if(ioflag.gt.0) then
        io=0
      else
        io=51
      end if
c
      open(48, file='stop.in', status='old')
      read(48,*) iistop
      if(iistop.gt.0) stop
c
      write(io, fmt="(/'Enter input file name')")
      read(5,85) file1
85    format(a72)
      write(io, fmt="(/'Enter output file name')")
      read(5,85) file2
c
      open(10, file=file1, form='unformatted', status='old')
      open(35, file='for.header', status='old')
      open(36, file='inv.header', status='old')
c
c     read in model dimensions
c
      read(35,115) xmin,xmax,ymin,ymax,zmin,zmax,size,nx,ny,nz
115   format(7f10.3,3i10)
c
      if(nx.gt.nxmax.or.ny.gt.nymax.or.nz.gt.nzmax) then
        write(0,9)
9       format(/'***  model is too large  ***'/)
        stop
      end if
c 
      nodestotal=nx*ny*nz 
c
      read(36,*) nxi,nyi,nzi
c
      if(nxi.lt.1.or.nyi.lt.1.or.nzi.lt.1.or.
     +   nxi.gt.nximax.or.nyi.gt.nyimax.or.nzi.gt.nzimax) then
        write(0,8)
8       format(/'***  original file has invalid size  ***'/)
        stop
      end if
c
      cellstotal=nxi*nyi*nzi
c
      write(io,25) xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,nodestotal 
25    format(/'model dimensions:'/ 
     +        '-----------------'/ 
     +        'xmin, xmax: ',2f10.3/'ymin, ymax: ',2f10.3/, 
     +        'zmin, zmax: ',2f10.3/ 
     +        'number of model nodes in each direction (x,y,z):',3i4/ 
     +        'total number of model nodes:                    ',i12) 
      write(io,35) size 
35    format('node spacing: ',f7.3,' km') 
      write(io,45) nxi,nyi,nzi,cellstotal
45    format(
     +'original number of cells in each direction (x,y,z):  ',3i4/
     +'original total number of cells:                      ',i12)
c
      ixinc=(nx-1)/nxi
      iyinc=(ny-1)/nyi
      izinc=(nz-1)/nzi
      num=(ixinc+1)*(iyinc+1)*(izinc+1)
c
      do 20 k=1,nz
         do 20 j=1,ny
            read(10) (ivel(i,j,k),i=1,nx)
20    continue
c
      close(10)
      open(11, file=file2, form='unformatted')
c
      do 30 k=1,nzi
         kmin=(k-1)*izinc+1         
         kmax=kmin+izinc
         do 30 j=1,nyi
            jmin=(j-1)*iyinc+1         
            jmax=jmin+iyinc
            do 30 i=1,nxi
               imin=(i-1)*ixinc+1         
               imax=imin+ixinc
               sum=0.
               do 400 ii=imin,imax
                  do 400 jj=jmin,jmax
                     do 400 kk=kmin,kmax
400                     sum=sum+1000./ivel(ii,jj,kk)
               smod(i,j,k)=sum/num
30    continue
c
      do 40 k=1,nzi
         do 40 j=1,nyi
40          write(11) (smod(i,j,k),i=1,nxi)
c
      stop
      end
