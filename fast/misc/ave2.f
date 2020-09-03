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
c        10 -- input:  input 3D data grid
c
c        12 -- output: ave 3D data grid
c
c        35 -- input:  parameters describing size of 3D grid
c
c
c     ----------------------------------------------------------------
c
c
      include 'ray.par'
c
      character*72 file1
      integer*2 veli1(nxmax,nymax,nzmax),veli2(nxmax,nymax,nzmax)
c
      write(*, fmt="(/'Enter input file name')")
      read(5,85) file1
85    format(a72)
      write(*, fmt="(/'Enter min and max depth (km)')")
      read(5,*) zamin,zamax
c
      open(10, file=file1, form='unformatted', status='old')
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
      nz1=nint((zamin-zmin)/size)+1
      nz2=nint((zamax-zmin)/size)+1
      naz=nz2-nz1+1
      write(0,*) nz1,nz2,naz
c
      do 10 k=1,nz
         do 10 j=1,ny
10           read(10) (veli1(i,j,k),i=1,nx)
c
      do 20 j=1,ny
         do 20 i=1,nx
             sum=0.
             do 30 k=nz1,nz2
30              sum=sum+veli1(i,j,k)
             iave=nint(sum/naz)
             do 40 k=1,nz
40              veli2(i,j,k)=iave
20    continue
c
      do 50 k=1,nz
         do 50 j=1,ny
50           write(12) (veli2(i,j,k),i=1,nx)
c
      stop
      end
