c
c     version 1.0  Jul 1996
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |                *****  V E L I N T 4  ****                    |
c     |                                                              |
c     |  Create the starting model given the bathymetry and 1D model |
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
c        10 -- input:  input 1D velocity model
c
c        11 -- output:  output 3D velocity 
c
c        12 -- input:  bathymetry
c
c        35 -- input:  parameters describing size of 3D grid
c
c
c     ----------------------------------------------------------------
c
c
      include 'ray.par'
c
      real d1(nzmax),d2(nzmax),v1(nzmax),v2(nzmax)
      integer*2 vel(nxmax,nymax,nzmax)
      real bath(nxmax,nymax)
c
      open(10, file='velint.in', status='old')
      open(11, file='vel.int', form='unformatted')
      open(12, file='bathymetry', form='unformatted', status='old')
      open(35, file='for.header', status='old')
c
c     read in model dimensions
c
      read(35,115) xmin,xmax,ymin,ymax,zmin,zmax,size,nx,ny,nz
115   format(7f10.3,3i10)
c
      if(nx.gt.nxmax.or.ny.gt.nymax) then
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
      size2=size/2.
      size4=size/4.
c
      i=1
100   read(10,*,end=999) d1(i),v1(i)
      read(10,*) d2(i),v2(i)
      v1(i)=v1(i)*1000.
      v2(i)=v2(i)*1000.
      i=i+1
      go to 100
999   nl=i-1
c
      do 20 j=1,ny
20       read(12) (bath(i,j),i=1,nx)
c
      do 30 i=1,nx
         do 30 j=1,ny
            db=bath(i,j)/1000.
            d2(1)=db
            d1(2)=db
c
            do 40 k=1,nz
               depth=zmin+float(k-1)*size
               du=depth-size2
               if(du.lt.zmin) du=zmin
               dl=depth+size2
               if(dl.gt.zmax) dl=zmax
               du4=depth-size4
               if(du4.lt.zmin) du4=zmin
               dl4=depth+size4
               if(dl4.gt.zmax) dl4=zmax
               do 50 l=1,nl
                  if(depth.ge.d1(l).and.depth.lt.d2(l)) then
                    vdepth=(v2(l)-v1(l))/(d2(l)-d1(l))*
     +                         (depth-d1(l))+v1(l)
                    go to 41
                  end if
50             continue
               go to 49
41             do 51 l=1,nl
                  if(du.ge.d1(l).and.du.lt.d2(l)) then
                    vu=(v2(l)-v1(l))/(d2(l)-d1(l))*
     +                       (du-d1(l))+v1(l)
                    go to 42
                  end if
51             continue
               go to 49
42             do 52 l=1,nl
                  if(dl.ge.d1(l).and.dl.lt.d2(l)) then
                    vl=(v2(l)-v1(l))/(d2(l)-d1(l))*
     +                       (dl-d1(l))+v1(l)
                    go to 43
                  end if
52             continue
               go to 49
43             do 53 l=1,nl
                  if(du4.ge.d1(l).and.du4.lt.d2(l)) then
                    vu4=(v2(l)-v1(l))/(d2(l)-d1(l))*
     +                       (du4-d1(l))+v1(l)
                    go to 44
                  end if
53             continue
               go to 49
44             do 54 l=1,nl
                  if(dl4.ge.d1(l).and.dl4.lt.d2(l)) then
                    vl4=(v2(l)-v1(l))/(d2(l)-d1(l))*
     +                       (dl4-d1(l))+v1(l)
                    go to 48
                  end if
54             continue
49             write(0,*) '***  could not find layer  ***'
               write(0,*) i,j,k,depth
               stop
48             vel(i,j,k)=(vu+2.*vu4+2.*vdepth+2.*vl4+vl)/8.
40          continue
30    continue
c
      do 60 k=1,nz
         do 60 j=1,ny
60          write(11) (vel(i,j,k),i=1,nx)
c
      stop
      end
