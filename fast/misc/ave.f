c
c     version 1.0  May 1995
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |                     *****  A V E  ****                       |
c     |                                                              |
c     |            Calculate average of two 3D data grids            |
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
      real vel1(nxmax),vel2(nxmax)
      character*72 file1,file2
      character*1 r4
      integer*2 veli1(nxmax),veli2(nxmax)
c
      write(*, fmt="(/'Enter file name 1')")
      read(5,85) file1
85    format(a72)
      write(*, fmt="(/'Enter file name 2')")
      read(5,85) file2
      write(*, fmt="(/'Enter  r  for r*4 file (default is i*2)')")
      read(5,95) r4
95    format(a1)
c
      open(10, file=file1, form='unformatted', status='old')
      open(11, file=file2, form='unformatted', status='old')
      open(12, file='ave.out', form='unformatted')
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
      nk=0
c 
      do 10 k=1,nz
         do 11 j=1,ny
           if(r4.eq.'r') then
             read(10,end=999) (vel1(i),i=1,nx)
             read(11,end=999) (vel2(i),i=1,nx)
             write(12) ((vel1(i)+vel2(i))/2.,i=1,nx)
           else
             read(10,end=999) (veli1(i),i=1,nx)
             read(11,end=999) (veli2(i),i=1,nx)
             do 20 i=1,nx
20              veli1(i)=nint((veli1(i)+veli2(i))/2.)
             write(12) (veli1(i),i=1,nx)
           end if
11       continue
         nk=nk+1
10    continue
c
      stop
c
999   if(nk.eq.1) then
        write(6,155) nk
155     format(/'>>> data file contains ',i3,' surface'/)
      else
        write(6,156) nk
156     format(/'>>> data file contains ',i3,' surfaces'/)
      end if
c
      stop
      end
