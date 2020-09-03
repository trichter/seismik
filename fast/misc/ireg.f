c
c     Determine the top cells at which regularization can
c     begin below an interface
c
      include 'ray.par'
c
      real boundary(nxmax,nymax)
      integer ireg(nximax,nyimax)
c
      open(35, file='for.header', status='old')
      open(36, file='inv.header', status='old')
      open(11, file='bathymetry', form='unformatted', status='old')
      open(12, file='ireg.cell', form='unformatted')
c
      write(6,335)
335   format('IREG: resample interface on inverse grid')
      open(50, file='nowrite', status='old', err=9999, iostat=ioflag)
9999  if(ioflag.gt.0) then
        io=0
      else
        io=51
      end if
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
     +   nxi.ge.nxmax.or.nyi.ge.nymax.or.nzi.ge.nzmax) then
        write(0,8)
8       format(/'***  original file has invalid size  ***'/)
        stop
      end if
c
      nodestotali=nxi*nyi*nzi
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
      write(io,45) nxi,nyi,nzi,nodestotali
45    format(
     +'original number of cells in each direction (x,y,z):  ',3i4/
     +'original total number of cells:                      ',i12)
c
      zinc=(zmax-zmin)/float(nzi)
c
      do j=1,ny
         read(11) (boundary(i,j),i=1,nx)
         do 10 i=1,nx
10          boundary(i,j)=boundary(i,j)/1000.
      enddo    
c
      ixf=(nx-1)/nxi
      iyf=(ny-1)/nyi
      do 20 i=1,nxi
         i1=(i-1)*ixf+1
         i2=i*ixf+1
         do 20 j=1,nyi
            j1=(j-1)*iyf+1
            j2=j*iyf+1
            zimax=zmin
            do 30 ii=i1,i2
               do 30 jj=j1,j2
                  if(boundary(ii,jj).gt.zimax) zimax=boundary(ii,jj)
30          continue
         ireg(i,j)=int((zimax-zmin)/zinc)+1    
20    continue
c
      do j=1,nyi
         write(12) (ireg(i,j),i=1,nxi)
      enddo
c
      stop
      end
