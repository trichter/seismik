c
c     version 1.0  Jul 1994
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |            **********  U P D A T E  **********               |   
c     |                                                              |
c     |            Update the 3D velocity model given                |
c     |                the slowness perturbation                     |
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
c        17 -- input:  slowness perturbation 
c
c        35 -- input:  parameters describing size of 3D grid
c
c
c     ----------------------------------------------------------------
c
c
      include 'ray.par'
c
      real pert(nxmax)
      integer*2 vel(nxmax),velnew(nxmax)
c
      write(6,335)
335   format('UPDATE: update velocity model')
      open(50, file='nowrite', status='old', err=9999, iostat=ioflag)
9999  if(ioflag.gt.0) then
        io=0
      else
        io=51
      end if
c
      open(48, file='stop.in', status='old')
      read(48,*) iistop
      if(iistop.gt.0) stop
c
      open(10, file='vel.mod', form='unformatted', 
     +     status='old')
      open(11, file='vel.new', form='unformatted')
      open(17, file='sl.pert', form='unformatted', 
     +     status='old')
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
      do 10 k=1,nz
         do 10 j=1,ny
           read(17) (pert(i),i=1,nx)
           read(10) (vel(i),i=1,nx)
           do 20 i=1,nx
20            velnew(i)=nint(1000./(1000./float(vel(i))+pert(i)))
           write(11) (velnew(i),i=1,nx)
10    continue
c
      stop
      end
