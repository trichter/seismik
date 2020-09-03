c
c     version 1.0  Jul 1998
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |                  *****  B O U N D  ****                      |
c     |                                                              |
c     |      Determine which velocities have hit their bounds        |
c     |                                                              |
c     |                 Written by Colin A. Zelt                     |
c     |                                                              |
c     |                Dept of Geology & Geophysics                  |
c     |                      Rice University                         |
c     |                  Houston, TX  77005-1892                     |
c     |                                                              |
c     ----------------------------------------------------------------
c
c
c     I/O units:
c
c        10 -- input:  input 3D data grid
c
c        11 -- output:  output 3D data grid
c
c        12 -- input: bathymetry
c
c        35 -- input:  parameters describing size of 3D grid
c
c
c     ----------------------------------------------------------------
c
c
      include 'ray.par'
c
      parameter(npair=10)
      character*72 filei,fileo,fileb
      integer*2 vel(nxmax,nymax,nzmax),velb(nxmax,nymax,nzmax)
      real bath(nxmax,nymax),
     +     v0(npair)/npair*-1./,vmin(npair),vmax(npair)
c
      namelist /rspar/ interface,v0,vmin,vmax,filei,fileo,fileb
c
      interface=0
c
      open(14, file='bound.in',status='old')
c
      read(14,rspar)
c
      open(10, file=filei, form='unformatted',status='old')
      open(11, file=fileo, form='unformatted')
      open(13, file=fileb, form='unformatted',status='old')
      open(35, file='for.header', status='old')
c
c     read in model dimensions
c
      read(35,115) xmin,xmax,ymin,ymax,zmin,zmax,size,nx,ny,nz
115   format(7f10.3,3i10)
c
      if(nx.gt.nxmax.or.ny.gt.nymax) then
        write(0,9)
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
      if(interface.eq.1) then
        open(12, file='bathymetry', form='unformatted',status='old')
        do 10 j=1,ny
10         read(12) (bath(i,j),i=1,nx)
      end if
c
      nv=0
      do i=1,npair
         if(v0(i).lt.0.) go to 1998
         v0(i)=1000.*v0(i)
         vmin(i)=1000.*vmin(i)
         vmax(i)=1000.*vmax(i)
         nv=nv+1
      enddo
1998  if(nv.eq.0) then
        write(0,*) '***  must specify v0, vmin and vmax values  ***'
        stop
      end if
c
      nfix=0
c
      do 20 k=1,nz
         do 20 j=1,ny
            read(10) (vel(i,j,k),i=1,nx)
20          read(13) (velb(i,j,k),i=1,nx)
c
      do 30 i=1,nx
         do 30 j=1,ny 
            if(interface.eq.1) then
              ns1=nint((bath(i,j)/1000.-zmin)/size)+1
              do k=1,ns1
                  vel(i,j,k)=0
              enddo
            else
              ns1=0
            end if
c
            if(ns1.lt.nz) then
              do 32 k=ns1+1,nz
c
               iflag=0
               do l=2,nv
                  if(float(velb(i,j,k)).le.v0(l).or.l.eq.nv) then
                    vpmin=(vmin(l)-vmin(l-1))/(v0(l)-v0(l-1))*
     +                    (velb(i,j,k)-v0(l-1))+vmin(l-1)
                    vpmax=(vmax(l)-vmax(l-1))/(v0(l)-v0(l-1))*
     +                    (velb(i,j,k)-v0(l-1))+vmax(l-1)
                    if(float(vel(i,j,k)).le.vpmin) then
                      iflag=1
                      go to 1999
                    end if
                    if(float(vel(i,j,k)).ge.vpmax) then
                      iflag=1
                    end if
                    go to 1999
                  end if
               enddo
c
1999           if(iflag.eq.1) then
                 nfix=nfix+1
                 vel(i,j,k)=1
               else
                 vel(i,j,k)=0
               end if

32            continue
            end if
30    continue
c
      do 40 k=1,nz
         do 40 j=1,ny
40          write(11) (vel(i,j,k),i=1,nx)
c
      write(6,45) nfix
45    format('number of velocities set to bounds: ',i10)
c
      stop
      end
