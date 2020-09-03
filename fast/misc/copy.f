c
c     version 1.0  Jan 1994
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |                    *****  C O P Y  ****                      |
c     |                                                              |
c     |      Archive and copy new velocity model to current model    |
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
c        10 -- input:  new velocity model 
c
c        11 -- output: current velocity model
c
c        12 -- output: archived velocity model
c
c        13 -- input:  archived velocity file names
c
c        14 -- input: current iteration number
c
c        35 -- input:  parameters describing size of 3D grid
c
c
c     ----------------------------------------------------------------
c
c
      include 'ray.par'
c
      integer*2 vel(nxmax)
      character filea*80,filen*10
c
      write(6,335)
335   format('COPY: archive and copy new velocity model')
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
      filen='velfile.  '
c
      open(10, file='vel.new', form='unformatted', status='old')
      open(11, file='vel.mod', form='unformatted', status='old')
      open(13, file='file.names', status='old')
      open(14, file='current.iteration', status='old')
      open(35, file='for.header', status='old')
c
      read(14,*) iteration
      rewind(14)
      write(14,*) iteration+1
c
      do 20 i=1,iteration
         read(13,5,end=999) filea
5        format(a80)
20    continue
c
      open(12, file=filea, form='unformatted')
c
      go to 98
c
999   write(io,85)
85    format(/'***  premature end of file in file.names  ***'/)
      id1=iteration/10
      id2=iteration-id1*10
      filen(9:10)=char(id1+48)//char(id2+48)
c
      open(12, file=filen, form='unformatted')
c
c     read in model dimensions
c
98    read(35,115) xmin,xmax,ymin,ymax,zmin,zmax,size,nx,ny,nz
115   format(7f10.3,3i10)
c
      if(nx.gt.nxmax) then
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
           read(10) (vel(i),i=1,nx)
           write(11) (vel(i),i=1,nx)
10         write(12) (vel(i),i=1,nx)
c
      stop
c
      end
