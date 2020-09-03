c
c     version 1.0  Nov 1993
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |            ********  A N O M A L Y  **********               |   
c     |                                                              |
c     |        Add a relative sinusoidal velocity anomaly            |   
c     |                  to a 3-D velocity model                     |
c     |                                                              |
c     |                 Written by Colin A. Zelt                     |   
c     |                                                              |
c     |                   Bullard Laboratories                       |
c     |                  University of Cambridge                     |   
c     |                  Cambridge, UK  CB3 0EZ                      |
c     |                                                              |
c     ----------------------------------------------------------------
c
c
c     I/O units:
c
c        10 -- input:  current velocity model
c
c        11 -- output: updated velocity model
c
c        35 -- input:  parameters describing size of 3-D grid
c
c
c     ----------------------------------------------------------------
c
c
      include 'ray.par'
c
      character*80 file0
      integer*2 vel(nxmax,nymax,nzmax)
c
      open(unit=11, file='vel.anomaly', form='unformatted')
      open(unit=35, file='for.header', status='old')
c
c     read in model dimensions
c
      read(35,115) xmin,xmax,ymin,ymax,zmin,zmax,size,nx,ny,nz
115   format(7f10.3,3i10)
c
      if(nx.gt.nxmax.or.ny.gt.nymax.or.nz.gt.nzmax) then
        write(6,9)
9       format(/'***  model is too large  ***'/)
        stop
      end if
c 
      nodestotal=nx*ny*nz 
      write(6,25) xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,nodestotal 
25    format(/'model dimensions:'/ 
     +        '-----------------'/ 
     +        'xmin, xmax: ',2f10.3/'ymin, ymax: ',2f10.3/, 
     +        'zmin, zmax: ',2f10.3/ 
     +        'number of model nodes in each direction (x,y,z):',3i4/ 
     +        'total number of model nodes:                    ',i12) 
      write(6,35) size 
35    format('node spacing: ',f7.3,' km') 
c
      write(6,*) 'Enter background velocity model file'
      read(5,*) file0
      write(6,*) 'Enter x and y anomaly width (km)'
      read(5,*) xsize,ysize
      write(6,*) 'Enter xshift and yshift (km)'
      read(5,*) xshift,yshift
      write(6,*) 'Enter theta (degrees)'
      read(5,*) theta
      theta=theta/57.29578
      sint=sin(theta)
      cost=cos(theta)
      write(6,*) 'Enter polarity'
      read(5,*) ipol
      write(6,*) 'Enter zmin and zmax of anomaly (km)'
      read(5,*) zamin,zamax
      write(6,*) 'Enter relative velocity anomaly (%)'
      read(5,*) vrel
      if(ipol.ne.-1) ipol=1
      vrel=ipol*vrel/100.
c
      if(zamin.lt.zmin) zamin=zmin
      if(zamax.gt.zmax) zamax=zmax
      nz1=nint((zamin-zmin)/size)+1
      nz2=nint((zamax-zmin)/size)+1
      xf=3.1415927/xsize
      yf=3.1415927/ysize
c
      open(unit=10, file=file0, form='unformatted', 
     +     status='old')
c
      do 10 k=1,nz
         do 10 j=1,ny
10         read(10) (vel(i,j,k),i=1,nx)
c
      do 40 j=1,ny
         y=ymin+float(j-1)*size
         do 40 i=1,nx
           x=xmin+float(i-1)*size
           xr=x*cost+y*sint
           yr=y*cost-x*sint
           an=vrel*sin(xf*(xr-xshift))*sin(yf*(yr-yshift))
           pan=1.+an
           do 20 k=nz1,nz2
              vel(i,j,k)=nint(float(vel(i,j,k))*pan)
20         continue
40    continue
c
999   do 50 k=1,nz
         do 50 j=1,ny
50         write(11) (vel(i,j,k),i=1,nx)
c
      stop
      end
