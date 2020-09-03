c
c     version 1.0  May 1995
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |                     *****  A V E  ****                       |
c     |                                                              |
c     |            Calculate average of four 3D data grids           |
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
c        patched for g77 compilation
c
c        Drew Brenders (brenders@geoladm.geol.queensu.ca)
c
c     I/O units:
c
c        10 -- input:  first 3D data grid
c
c        11 -- input:  second 3D data grid
c
c        12 -- output: ave (file1 + file2)/2
c
c        35 -- input:  parameters describing size of 3D grid
c
c
c     ----------------------------------------------------------------
c
c
      include 'ray.par'
c
      character*72 file1,file2,file3,file4
      integer*2 vel(4,nxmax),vela(nxmax)
c
      write(*, fmt="(/'Enter file name 1')")
      read(5,85) file1
85    format(a72)
      write(*, fmt="(/'Enter file name 2')")
      read(5,85) file2
      write(*, fmt="(/'Enter file name 3')")
      read(5,85) file3
      write(*, fmt="(/'Enter file name 4')")
      read(5,85) file4
c
      open(11, file=file1, form='unformatted', status='old')
      open(12, file=file2, form='unformatted', status='old')
      open(13, file=file3, form='unformatted', status='old')
      open(14, file=file4, form='unformatted', status='old')
      open(15, file='ave.out', form='unformatted')
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
      do 10 k=1,nz
         do 11 j=1,ny
           read(11) (vel(1,i),i=1,nx)
           read(12) (vel(2,i),i=1,nx)
           read(13) (vel(3,i),i=1,nx)
           read(14) (vel(4,i),i=1,nx)
           do 20 i=1,nx
20            vela(i)=nint((vel(1,i)+vel(2,i)+vel(3,i)+vel(4,i))/4.)
           write(15) (vela(i),i=1,nx)
11       continue
10    continue
c
      stop
      end
