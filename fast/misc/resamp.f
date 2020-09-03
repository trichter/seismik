c
c     version 1.0  Mar 1995
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |              **********  R E S A M P  *********              |
c     |                                                              |
c     |        resample a 3D file onto a different grid spacing      |
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
c        patched for g77 compilation
c
c        Drew Brenders (brenders@geoladm.geol.queensu.ca)
c
c     I/O units:
c
c        10 -- input:  3D integer*2 file 
c
c        11 -- output: 3D integer*2 file sampled on new grid
c
c        35 -- input:  parameters describing size of 3D grid
c
c
c     ----------------------------------------------------------------
c
c
      include 'ray.par'
c
      integer*2 vel(nxmax,nymax,nzmax),datai(nxmax)
      character filei*72,fileo*72
c
      write(*, fmt="(/'Enter input file name')")
      read(5,85) filei
85    format(a72)
      write(*, fmt="(/'Enter output file name')")
      read(5,85) fileo
c
      open(10, file=filei, form='unformatted', status='old')
      open(11, file=fileo, form='unformatted')
      open(35, file='for.header', status='old')
c
c     read in model dimensions
c
98    read(35,115) xmin,xmax,ymin,ymax,zmin,zmax,size,nx,ny,nz
115   format(7f10.3,3i10)
c
      if(nx.gt.nxmax) then
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
      write(6,1)
1     format(/'enter new grid spacing')
      read(5,*) size2
c
      do 10 k=1,nz
         do 10 j=1,ny
           read(10) (datai(i),i=1,nx)
           do 20 i=1,nx
20            vel(i,j,k)=datai(i)
10    continue      
c
      nx2=int(1.0001*(xmax-xmin)/size2)+1
      ny2=int(1.0001*(ymax-ymin)/size2)+1
      nz2=int(1.0001*(zmax-zmin)/size2)+1
      xmax2=xmin+float(nx2-1)*size2
      ymax2=ymin+float(ny2-1)*size2
      zmax2=zmin+float(nz2-1)*size2
c
      write(0,125) xmin,xmax2,ymin,ymax2,zmin,zmax2,size2,
     +             nx2,ny2,nz2
125   format(/'new model header line:'/7f10.3,3i10)
c
      do 30 k=1,nz2
         zp=zmin+float(k-1)*size2
         do 30 j=1,ny2
            yp=ymin+float(j-1)*size2
            do 40 i=1,nx2
               xp=xmin+float(i-1)*size2
               ixc=int((xp-xmin)/size)+1
               iyc=int((yp-ymin)/size)+1
               izc=int((zp-zmin)/size)+1
c
               v1=vel(ixc,iyc,izc)
               v2=vel(ixc+1,iyc,izc)
               v3=vel(ixc,iyc+1,izc)
               v4=vel(ixc+1,iyc+1,izc)
               v5=vel(ixc,iyc,izc+1)
               v6=vel(ixc+1,iyc,izc+1)
               v7=vel(ixc,iyc+1,izc+1)
               v8=vel(ixc+1,iyc+1,izc+1)
c
               xfrac=xp-(float(ixc-1)*size+xmin)
               yfrac=yp-(float(iyc-1)*size+ymin)
               zfrac=zp-(float(izc-1)*size+zmin)
c
               if(ny.eq.1) then
                 v3=v1
                 v4=v2
                 v7=v5
                 v8=v6
                 yfrac=0.
               end if
c
               vleft= (size*(v3-v1)*yfrac+size*(v5-v1)*zfrac+
     +                (v1-v3+v7-v5)*yfrac*zfrac)/size2+v1
               vright=(size*(v4-v2)*yfrac+size*(v6-v2)*zfrac+
     +                (v2-v4+v8-v6)*yfrac*zfrac)/size2+v2
c
               datai(i)=((size-xfrac)*vleft+xfrac*vright)/size
40          continue
c
            write(11) (datai(i),i=1,nx2)
30    continue
c
      stop
c
      end
