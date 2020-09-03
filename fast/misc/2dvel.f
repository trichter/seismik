c
c     version 1.0  Jun 1997
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |                  *****  2 D V E L  ****                      |
c     |                                                              |
c     |        calculate the average 2D model from a 3D model        |
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
c        10 -- input:  input 3D model
c
c        11 -- output:  output 3D model (that is the 2D average model)
c
c        12 -- input: num.cell   
c
c        13 -- output:  2D model ascii file
c
c        35 -- input:  parameters describing size of 3D grid
c
c
c     ----------------------------------------------------------------
c
c
      include 'ray.par'
c
      character*72 filei,fileo
      integer*2 vel(nxmax,nymax,nzmax),velo(nxmax,nymax,nzmax),
     +          vave
c
      write(*, fmt="(/'Enter input file name')")
      read(5,85) filei
85    format(a72)
      write(*, fmt="(/'Enter output file name (default is input)')")
      read(5,85) fileo
      if(fileo.eq.'') fileo=filei
c
      open(10, file=filei, form='unformatted', status='old')
      open(11, file=fileo, form='unformatted')
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
      fny=float(ny)
c
      do 30 k=1,nz
         do 30 i=1,nx
            sum=0.
            sum2=0.
            do 40 j=1,ny 
40             sum=sum+vel(i,j,k) 
         fvave=sum/fny
         vave=nint(sum/fny)
         do 50 j=1,ny
50          velo(i,j,k)=vave
30    continue
c
      do 60 k=1,nz
         do 60 j=1,ny
60          write(11) (velo(i,j,k),i=1,nx)
c
      stop
      end
