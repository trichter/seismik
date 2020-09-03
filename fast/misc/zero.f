c
c     version 1.0  Apr 1995
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |               ********  Z E R O  **********                  |
c     |                                                              |
c     |           Create the files necessary for slowness            |
c     |          tomography or calculating spread function           |
c     |               and initialize with zero values                |
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
c        17 -- output: either file sl.sums or res.ker
c
c        18 -- output: file weightcell.ray
c
c        19 -- output: file numcell.ray
c
c        20 -- output: file depth.sums
c
c        35 -- input:  parameters describing size of 3D grid
c
c
c     ----------------------------------------------------------------
c
c
      include 'ray.par'
c
      real rarr(nxmax)
      integer*2 iarr(nxmax)
c
      write(6,335)
335   format('ZERO: create and initialize sum and count files')
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
      write(io,1)
1     format(/
     +'enter  1  for slowness tomography'/
     +'   or  2  for spread function (resolution kernel)')
      read(5,*) itype
c
      if(itype.eq.1) then
        open(17, file='sl.sums', form='unformatted')
        open(18, file='weight.cell', form='unformatted')
        open(19, file='num.cell', form='unformatted')
      end if
      if(itype.eq.2) 
     +  open(17, file='res.ker', form='unformatted')
c
      open(35, file='inv.header', status='old')
c
c     read in model dimensions
c
      read(35,*) nx,ny,nz
c 
      nodestotal=nx*ny*nz 
      write(io,25) nx,ny,nz,nodestotal 
25    format(/'number of model nodes in each direction (x,y,z):',3i4/ 
     +        'total number of model nodes:                    ',i12) 
c
      if(nx.gt.nxmax.or.ny.gt.nymax.or.nz.gt.nzmax) then
        write(0,3)
3       format(/'***  model is too large  ***'/)
        stop
      end if

c
      do 10 i=1,nx
         iarr(i)=0
10       rarr(i)=0.
c
      if(itype.eq.1) then
        do 20 k=1,nz
           do 20 j=1,ny
              write(17) (rarr(i),i=1,nx)
              write(18) (rarr(i),i=1,nx)
20            write(19) (iarr(i),i=1,nx)
      end if
c
      if(itype.eq.2) then
        do 50 k=1,nz
           do 50 j=1,ny
50            write(17) (rarr(i),i=1,nx)
      end if
c
      stop
      end
