c
c     version 1.0  Nov 1993
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |            ********  A N O M A L Y  **********               |   
c     |                                                              |
c     |            Add a rectangular velocity anomaly                |   
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
      integer*2 vel(nxmax,nymax,nzmax),vhold(nxmax,nymax,nzmax)
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
     +'enter xmin,xmax,ymin,ymax,zmin,zmax of velocity anomaly')
      read(5,*) xamin,xamax,yamin,yamax,zamin,zamax
      if(xamin.eq.xamax.and.yamin.eq.yamax.and.zamin.eq.zamax) go to 999
      write(6,2)
2     format(/'enter velocity anomaly to add and 0 for absolute,'/
     +        '1 for relative additive and 2 for relative percentage'/
     +        ' (type=3 to stop)')
      read(5,*) vanomaly,itype
c
      if(itype.eq.3) go to 999
c
      write(6,3)
3     format(/'enter thickness of transition zone in nodes')
      read(5,*) ntrans
c
      if(xamin.lt.xmin) xamin=xmin
      if(xamax.gt.xmax) xamax=xmax
      if(yamin.lt.ymin) yamin=ymin
      if(yamax.gt.ymax) yamax=ymax
      if(zamin.lt.zmin) zamin=zmin
      if(zamax.gt.zmax) zamax=zmax
c
      nx1=nint((xamin-xmin)/size)+1
      nx2=nint((xamax-xmin)/size)+1
      ny1=nint((yamin-ymin)/size)+1
      ny2=nint((yamax-ymin)/size)+1
      nz1=nint((zamin-zmin)/size)+1
      nz2=nint((zamax-zmin)/size)+1
c
      if(nx2-nx1.lt.2*ntrans.or.ny2-ny1.lt.2*ntrans.or.
     +   nz2-nz1.lt.2*ntrans) then
         write(0,*) nx1,nx2,ny1,ny2,nz1,nz2,2*ntrans
         write(6,5)
5        format(/
     +   '***  transition zone is too large for anomaly size  ***'/)
         go to 1000
      end if
c
      npts=0
c
      do 200 i=nx1,nx2
         do 200 j=ny1,ny2
            do 200 k=nz1,nz2
               npts=npts+1
200            vhold(i,j,k)=vel(i,j,k)
c
      if(itype.eq.0) then
        vanomaly=vanomaly*1000.
        if(ntrans.eq.0) then
          do 20 i=nx1,nx2
             do 20 j=ny1,ny2
                do 20 k=nz1,nz2
20                 vel(i,j,k)=vanomaly
        else
          do 21 l=1,ntrans+1
             do 22 i=nx1+l-1,nx2-l+1
                do 22 j=ny1+l-1,ny2-l+1
                   do 22 k=nz1+l-1,nz2-l+1
22                    vel(i,j,k)=(float(ntrans-l+1)*vhold(i,j,k)+
     +                         float(l)*vanomaly)/float(ntrans+1)
21        continue
        end if
      end if
c
      if(itype.eq.1) then
        vanomaly=vanomaly*1000.
        if(ntrans.eq.0) then
          do 30 i=nx1,nx2
             do 30 j=ny1,ny2
                do 30 k=nz1,nz2
30                 vel(i,j,k)=vel(i,j,k)+vanomaly
        else
          do 31 l=1,ntrans+1
             do 32 i=nx1+l-1,nx2-l+1
                do 32 j=ny1+l-1,ny2-l+1
                   do 32 k=nz1+l-1,nz2-l+1
32                    vel(i,j,k)=vel(i,j,k)+vanomaly/float(ntrans+1)
31        continue
        end if
      end if
c
      if(itype.eq.2) then
        if(ntrans.eq.0) then
          vanomaly=1.+vanomaly/100.
          do 40 i=nx1,nx2
             do 40 j=ny1,ny2
                do 40 k=nz1,nz2
40                 vel(i,j,k)=vel(i,j,k)*vanomaly
        else
          vanomaly=vanomaly/((ntrans+1)*100.)
          do 41 l=1,ntrans+1
             do 42 i=nx1+l-1,nx2-l+1
                do 42 j=ny1+l-1,ny2-l+1
                   do 42 k=nz1+l-1,nz2-l+1
42                    vel(i,j,k)=vel(i,j,k)+vhold(i,j,k)*vanomaly
41        continue
        end if
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
