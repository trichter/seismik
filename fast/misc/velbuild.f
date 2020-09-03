c
c     version 1.0  Mar 1995
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |               ******  V E L B U I L D  *****                 |   
c     |                                                              |
c     |        build a 1D, 2D or 3D velocity model for input to FD   |
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
c        10 -- input:  1D velocity-depth profile
c
c        12 -- output: 1D velocity-depth profile sampled on grid
c
c     20-22 -- input:  2D model from RAYINVR sampled on uniform grid
c
c        34 -- input:  parameters describing size of 3-D grid
c
c        36 -- output: 3D integer*2 velocity model
c
c
c     ----------------------------------------------------------------
c
c
      include 'ray.par'
c
      real velx1(nxmax),velx2(nxmax),zin(1000),vin(1000),
     +     vel(0:nxmax+1,0:nymax+1,0:nzmax+1)
      integer*2 veli(nxmax)
c
      write(6,2)
2     format(/'enter  1  to build from 1 RAYINVR model'/
     +        'enter  2  to build from 2 RAYINVR models'/
     +        'enter  3  to build from a 1D velocity-depth profile')
      read(5,*) imod
c
      open(unit=34, file='for.header', status='old')
c
c     read in model dimensions
c
      read(34,115) xmin,xmax,ymin,ymax,zmin,zmax,size,nx,ny,nz
115   format(7f10.3,3i10)
c
      if(nx.gt.nxmax.or.nz.gt.nzmax.or.ny.gt.nymax) then
        write(6,1)
1       format(/'***  model is too large  ***'/)
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
      if(imod.eq.1) then
c
c       read in single model from output of RAYINVR and expand it 
c       to 3D in y direction
c
        open(unit=20, file='v.in', status='old')
c
        read(20,5) xmin,xmax,zmin,zmax,size,nx,nz
5       format(5f10.3,2i10)
c     
        if(nx.gt.nxmax.or.nz.gt.nzmax) then
          write(6,1) 
          stop
        end if
c
        do 10 k=1,nz
           read(20,15) (velx1(i),i=1,nx)
15         format(10f10.3)
           do 11 i=1,nx
11            vel(i,1,k)=velx1(i)
10      continue
c
        ny=nint((ymax-ymin)/size)+1
c
        ymaxold=ymax
        ymax=ymin+float(ny-1)*size
c     
        if(ny.gt.nymax) then
          write(6,1) 
          stop
        end if
c
        do 20 i=1,nx
           do 20 k=1,nz
              do 20 j=2,ny
20               vel(i,j,k)=vel(i,1,k) 
c
      end if
c
      if(imod.eq.2) then
c
c       read in two models from output of RAYINVR and place them at
c       y=ymin and y=ymax and linearly interpolate between them
c
        open(unit=21, file='v1.in', status='old')
        open(unit=22, file='v2.in', status='old')
c
        read(21,5) xmin,xmax,zmin,zmax,size,nx,nz
        read(22,5) x2min,x2max,z2min,z2max,size2,nx2,nz2
        if(xmin.ne.x2min.or.xmax.ne.x2max.or.zmin.ne.z2min.or.
     +    zmax.ne.z2max.or.size.ne.size2.or.nx.ne.nx2.or.nz.ne.nz2) 
     +    then
          write(6,785)
785       format(/
     +    '***  imod=2 and two models are of different size  ***')
          stop
        end if
c     
        if(nx.gt.nxmax.or.nz.gt.nzmax) then
          write(6,1) 
          stop
        end if
c
        ny=nint((ymax-ymin)/size)+1
c
        ymaxold=ymax
        ymax=ymin+float(ny-1)*size
c
        if(ny.gt.nymax) then
          write(6,1)
          stop
        end if
c
        do 30 k=1,nz
           read(21,15) (velx1(i),i=1,nx)
           read(22,15) (velx2(i),i=1,nx)
           do 31 i=1,nx
              vel(i,1,k)=velx1(i)
31            vel(i,ny,k)=velx2(i)
30      continue
c
        do 40 i=1,nx
           do 40 k=1,nz
              do 40 j=2,ny-1
40               vel(i,j,k)=(float(ny-j)*vel(i,1,k)+float(j-1)*
     +                      vel(i,ny,k))/float(ny-1) 
c
      end if
c
      if(imod.eq.3) then
c
        open(10, file='vel1d.in', status='old')
        open(12, file='vel1d.out')
c
        i=1
100     read(10,*,end=999) zin(i),vin(i)
        i=i+1
        go to 100
999     npts=i-1
c
        do 200 k=1,nz
           depth=zmin+float(k-1)*size
           do 210 j=2,npts
              if(depth.ge.zin(j-1).and.depth.le.zin(j)) then
                velocity=(depth-zin(j-1))*(vin(j)-vin(j-1))/
     +                   (zin(j)-zin(j-1))+vin(j-1)
                go to 220
              end if
210        continue
           write(6,3)
3          format(/'***  1D profile not shallow or deep enough  ***'/)
           stop
220        do 230 i=1,nx
              do 230 j=1,ny
230              vel(i,j,k)=velocity
           write(6,4) depth,velocity
           write(12,4) depth,velocity
4          format(2f10.3)
200     continue 
c
      end if
c
      open(unit=36, file='vel.mod', form='unformatted')
      do 130 k=1,nz
         do 130 j=1,ny
            do 131 i=1,nx
131            veli(i)=vel(i,j,k)*1000.
130         write(36) (veli(i),i=1,nx)
c
      stop   
      end         
