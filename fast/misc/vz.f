c
c     version 1.0  Nov 1994
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |                  ********  V Z  *******                      |
c     |                                                              |
c     |      calculate 1D velocity-depth profile from a 3D model     |
c     |                                                              |
c     |                 Written by Colin A. Zelt                     |
c     |                                                              |
c     |                   Bullard Laboratories                       |
c     |                  University of Cambridge                     |
c     |                  Cambridge, UK  CB3 0EZ                      |
c     |                                                              |
c     ----------------------------------------------------------------
c
c        patched for g77 compilation
c
c        Drew Brenders (brenders@geoladm.geol.queensu.ca)
c
c     I/O units:
c
c        10 -- input:  input 3-D model
c
c        13 -- output:  1-D model ascii file
c
c        35 -- input:  parameters describing size of 3-D grid
c
c
c     ----------------------------------------------------------------
c
c
      include 'ray.par'
c
      character*72 filei
      integer*2 vel(nxmax,nymax,nzmax)
c
      write(*, fmt="(/'Enter input file name')")
      read(5,85) filei
85    format(a72)
      write(*, fmt="(/'Enter (x,y) position (km)')")
      read(5,*) xpos,ypos
c
      open(10, file=filei, form='unformatted', 
     +     status='old')
      open(13, file='vz.out')
      open(35, file='for.header', status='old')
c
c     read in model dimensions
c
      read(35,115) xmin,xmax,ymin,ymax,zmin,zmax,size,nx,ny,nz
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
      do 20 k=1,nz
         do 20 j=1,ny
20          read(10) (vel(i,j,k),i=1,nx)
c
      size2=size**2
      size1000=size*1000.
      zinc=(zmax-zmin)/float(nz-1)
      ixc=int((xpos-xmin)/size)+1
      iyc=int((ypos-ymin)/size)+1
      xfrac=xpos-(float(ixc-1)*size+xmin)
      yfrac=ypos-(float(iyc-1)*size+ymin)
c
      do 30 k=1,nz
         zpos=zmin+zinc*float(k-1)
         izc=k
c
         v1=vel(ixc,iyc,izc)
         v2=vel(ixc+1,iyc,izc)
         v3=vel(ixc,iyc+1,izc)
         v4=vel(ixc+1,iyc+1,izc)
         v5=vel(ixc,iyc,izc+1)
         v6=vel(ixc+1,iyc,izc+1)
         v7=vel(ixc,iyc+1,izc+1)
         v8=vel(ixc+1,iyc+1,izc+1)
         zfrac=zpos-(float(izc-1)*size+zmin)
c
         vleft= (size*(v3-v1)*yfrac+size*(v5-v1)*zfrac+
     +          (v1-v3+v7-v5)*yfrac*zfrac)/size2+v1
         vright=(size*(v4-v2)*yfrac+size*(v6-v2)*zfrac+
     +          (v2-v4+v8-v6)*yfrac*zfrac)/size2+v2
c
         velocity=((size-xfrac)*vleft+xfrac*vright)/size1000

         write(13,45) -zpos,velocity
45       format(2f10.3)
c
30    continue
c
      stop
      end
