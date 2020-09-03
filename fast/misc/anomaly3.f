c
c     version 1.0  Nov 1993
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |           ********  A N O M A L Y 3  **********              |   
c     |                                                              |
c     |            Add a disk-shaped velocity anomaly                |   
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
      integer*2 vel(nxmax,nymax,nzmax)
c
      open(unit=10, file='vel.mod', form='unformatted', 
     +     status='old')
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
      do 10 k=1,nz
         do 10 j=1,ny
10         read(10) (vel(i,j,k),i=1,nx)
c
1000  write(6,1)
1     format(/
     +'enter x0,y0,xradius,yradius,zmin,zmax of velocity anomaly')
      read(5,*) x0,y0,xradius,yradius,zamin,zamax
      write(6,2)
2     format(/'enter velocity anomaly to add and 0 for absolute,'/
     +        '1 for relative additive and 2 for relative percentage'/
     +        ' (velocity=0 to stop)')
      read(5,*) vanomaly,itype
c
      if(vanomaly.eq.0.) go to 999
c
      if(zamin.lt.zmin) zamin=zmin
      if(zamax.gt.zmax) zamax=zmax
c
      nz1=nint((zamin-zmin)/size)+1
      nz2=nint((zamax-zmin)/size)+1
c
      npts=0
c
      if(itype.eq.0) then
        vanomaly=vanomaly*1000.
        do 20 i=1,nx
           x=float(i-1)*size+xmin
           do 20 j=1,ny
              y=float(j-1)*size+ymin
              dist=((x-x0)/xradius)**2+((y-y0)/yradius)**2
              if(dist.le.1.) then
                do 21 k=nz1,nz2
                   npts=npts+1
21                 vel(i,j,k)=vanomaly
              end if
20         continue
      end if
c
      if(itype.eq.1) then
        vanomaly=vanomaly*1000.
        do 30 i=1,nx
           x=float(i-1)*size+xmin
           do 30 j=1,ny
              y=float(j-1)*size+ymin
              dist=((x-x0)/xradius)**2+((y-y0)/yradius)**2
              if(dist.le.1.) then
                do 31 k=nz1,nz2
                   npts=npts+1
31                 vel(i,j,k)=vel(i,j,k)+vanomaly
              end if
30         continue
      end if
c
      if(itype.eq.2) then
        vanomaly=1.+vanomaly/100.
        do 40 i=1,nx
           x=float(i-1)*size+xmin
           do 40 j=1,ny
              y=float(j-1)*size+ymin
              dist=((x-x0)/xradius)**2+((y-y0)/yradius)**2
              if(dist.le.1.) then
                do 41 k=nz1,nz2
                   npts=npts+1
41                 vel(i,j,k)=vel(i,j,k)*vanomaly
              end if
40         continue
      end if
c
      write(6,4) npts
4     format(/'number of points inside anomaly: ',i10/)
c
      go to 1000
c
999   do 50 k=1,nz
         do 50 j=1,ny
50         write(11) (vel(i,j,k),i=1,nx)
c
      stop
      end
