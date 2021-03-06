c
c     version 1.0  Mar 1995
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |               **********  H 2 Z  *********                   |
c     |                                                              |
c     |          convert a Hole real file to a Zelt integer*2 file   |
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
c        10 -- input:  real*4 file 
c
c        11 -- output: integer*2 file
c
c        35 -- input:  parameters describing size of 3D grid
c
c
c     ----------------------------------------------------------------
c
c
      include 'ray.par'
c
      integer*2 datai(nxmax)
      real data(nxmax)
      character filei*72,fileo*72
      dmin=1.e10
      dmax=-1.e10
      iflag=0
c
      write(*, fmt="(/'Enter input file name')")
      read(5,85) filei
85    format(a72)
      write(*, fmt="(/'Enter output file name')")
      read(5,85) fileo
c
      open(unit=11, file=fileo, form='unformatted')
      open(unit=35, file='for.header', status='old')
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
      open(unit=10, file=filei, form='unformatted', access='direct',
     +     recl=4*nx, status='old')
c
      write(6,1)
1     format(/'enter multiplier')
      read(5,*) xmult
c
      do 10 k=1,nz
         do 10 j=1,ny
           read(10,rec=((k-1)*ny+j)) (data(i),i=1,nx)
           do 20 i=1,nx
              if(data(i).lt.dmin) dmin=data(i)
              if(data(i).gt.dmax) dmax=data(i)
              di=data(i)*xmult
              if(abs(di).gt.32767.) iflag=1
20            datai(i)=nint(di)
           write(11) (datai(i),i=1,nx)
10    continue      
c
      write(6,5) dmin,dmax
5     format(/'dmin, dmax: ',2f15.5/)
c
      if(iflag.eq.1) write(6,55)
55    format(/'***  multiplier too large  ***'/)
c
      stop
c
      end
