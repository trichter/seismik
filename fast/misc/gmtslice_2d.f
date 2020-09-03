c
c     version 1.0  Apr 1996
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |             ********  G M T S L I C E  *****                 |
c     |                                                              |
c     |        Output slice of 2D volume for GMT plotting            |
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
c        10 -- input:  3D data grid
c
c        35 -- input:  parameters describing size of forward 3D grid
c
c
c     ----------------------------------------------------------------
c
c
      include 'ray.par'
c
      integer*2 vel(nxmax,nymax,nzmax),sample(nxmax,nymax,nzmax)
      character*72 file1
      character*9 file2
      real ninterp
c
      common /blk/ xmin,ymin,zmin,size,nx,ny,nz,vel,sample
c
      data file2/'y2d.slice'/
c
      write(0,335)
335   format('GMTSLICE: Output slices for GMT plotting')
      open(50, file='nowrite', status='old', err=9999, iostat=ioflag)
9999  if(ioflag.gt.0) then
        io=0
      else
        io=51
      end if
c
      write(io, fmt="(/'Enter input file name')")
      read(5,85) file1
85    format(a72)
      write(io, fmt="(/'Enter npts in x,z directions')")
      read(5,*) nxs,nzs
      write(io, fmt="(/'Enter data multiplier')")
      read(5,*) xmult
      write(io, fmt="(/'Enter  1  to mask unsampled cells')")
      read(5,*) mask
c
      open(10, file=file1, form='unformatted', status='old')
      if(mask.eq.1)
     +  open(12, file='num.cell', form='unformatted', status='old')
      open(35, file='for.header', status='old')
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
      write(io,25) xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,nodestotal 
25    format(/'model dimensions:'/ 
     +        '-----------------'/ 
     +        'xmin, xmax: ',2f10.3/'ymin, ymax: ',2f10.3/, 
     +        'zmin, zmax: ',2f10.3/ 
     +        'number of model nodes in each direction (x,y,z):',3i4/ 
     +        'total number of model nodes:                    ',i12) 
      write(io,35) size 
35    format('node spacing: ',f7.3,' km') 
c
      xinc=(xmax-xmin)/float(nxs-1)
      zinc=(zmax-zmin)/float(nzs-1)
c
      j=1
      do 20 k=1,nz
20          read(10) (vel(i,j,k),i=1,nx)
c
      if(mask.eq.1) then
        do 21 k=1,nz
21            read(12) (sample(i,j,k),i=1,nx)
      end if
      samp=1.
c
      close(10)
c
      open(11, file=file2)
      ypos=0.
      do 210 i=1,nxs
         xpos=xmin+float(i-1)*xinc
         do 210 k=1,nzs
            zpos=zmin+float(k-1)*zinc
            value=xmult*vinterp(xpos,ypos,zpos)
            if(mask.eq.1) samp=ninterp(xpos,ypos,zpos)
            if(samp.gt.0.) then
              write(11,*) xpos,zpos,value
            else
c             write(11,*) xpos,zpos,' NaN'
            end if
210      continue
      close(11)
c
      stop
      end
c
c     ----------------------------------------------------------------
c
      function vinterp(xp,yp,zp)
c
c     calculate the value at the point (xp,yp,zp), in general not
c     a grid point
c
      include 'ray.par'
c
      common /blk/ xmin,ymin,zmin,size,nx,ny,nz,vel,sample
c
      integer*2 vel(nxmax,nymax,nzmax),sample(nxmax,nymax,nzmax)
c
      ixc=min(nx,int((xp-xmin)/size)+1)
      iyc=1
      izc=min(nz,int((zp-zmin)/size)+1)
c
      v1=vel(ixc,iyc,izc)
      v2=vel(ixc+1,iyc,izc)
      v3=v1
      v4=v2
      v5=vel(ixc,iyc,izc+1)
      v6=vel(ixc+1,iyc,izc+1)
      v7=v5
      v8=v6
      xfrac=xp-(float(ixc-1)*size+xmin)
      yfrac=0.
      zfrac=zp-(float(izc-1)*size+zmin)

      vleft= (size*(v3-v1)*yfrac+size*(v5-v1)*zfrac+
     +       (v1-v3+v7-v5)*yfrac*zfrac)/(size*size)+v1
      vright=(size*(v4-v2)*yfrac+size*(v6-v2)*zfrac+
     +       (v2-v4+v8-v6)*yfrac*zfrac)/(size*size)+v2
c
      vinterp=((size-xfrac)*vleft+xfrac*vright)/size
c
      return
      end
c
c     ----------------------------------------------------------------
c
      function ninterp(xp,yp,zp)
c
c     calculate the value at the point (xp,yp,zp), in general not
c     a grid point
c
      include 'ray.par'
c
      common /blk/ xmin,ymin,zmin,size,nx,ny,nz,vel,sample
c
      real ninterp
      integer*2 vel(nxmax,nymax,nzmax),sample(nxmax,nymax,nzmax)
c
      ixc=min(nx,int((xp-xmin)/size)+1)
      iyc=1
      izc=min(nz,int((zp-zmin)/size)+1)
c
      v1=sample(ixc,iyc,izc)
      v2=sample(ixc+1,iyc,izc)
      v3=v1
      v4=v2
      v5=sample(ixc,iyc,izc+1)
      v6=sample(ixc+1,iyc,izc+1)
      v7=v5
      v8=v6
      xfrac=xp-(float(ixc-1)*size+xmin)
      yfrac=0.
      zfrac=zp-(float(izc-1)*size+zmin)

      vleft= (size*(v3-v1)*yfrac+size*(v5-v1)*zfrac+
     +       (v1-v3+v7-v5)*yfrac*zfrac)/(size*size)+v1
      vright=(size*(v4-v2)*yfrac+size*(v6-v2)*zfrac+
     +       (v2-v4+v8-v6)*yfrac*zfrac)/(size*size)+v2
c
      ninterp=((size-xfrac)*vleft+xfrac*vright)/size
c
      return
      end
