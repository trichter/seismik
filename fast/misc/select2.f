c
c     version 1.0  Jul 1996
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |            **********  S E L E C T  **********               |
c     |                                                              |
c     |                   Select the next model                      |
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
c        11 -- input:  chi.values
c
c        12 -- input:  lambda.values
c
c        13 -- input:  best velocity model 
c
c        14 -- output: vel.mod
c
c        15 -- output:  lambda value
c
c        16 -- input: input parameters
c
c        17 -- input: input parameters
c
c        35 -- input:  parameters describing size of 3D grid
c
c
c     ----------------------------------------------------------------
c
      include 'ray.par'
c
      integer*2 vel(nxmax)
      character vfile*6
      real lambda,lambmin
      integer chione
c
      namelist /selpar/ strfac
c
      strfac=1.
      chione=0
c
      vfile='vel.  '
      imin=0
      chidif=1.e20
c
      write(6,335)
335   format('SELECT: select the new model')
      open(50, file='nowrite', status='old', err=9999, iostat=ioflag)
9999  if(ioflag.gt.0) then
        io=0
      else
        io=51
      end if
c
      open(11, file='chi.values', status='old')
      open(12, file='lambda.values', status='old')
c
c      read(11,1) 
c      read(11,1) 
c      read(11,1) 
1     format(' ')
      read(11,*) 
      read(11,*) 
c
      i=1
100   read(11,5,end=99) chi
5     format(43x,f13.4)
c
      if(chi.lt.1.1) chione=1
c
      read(12,15,end=99) lambda
15    format(17x,e15.6)
      if(abs(chi-1.).lt.chidif) then
        chidif=abs(chi-1.)
        chimin=chi
        lambmin=lambda
        imin=i
      end if
c
      i=i+1
      go to 100
c
99    id1=imin/10
      id2=imin-id1*10
      if(imin.lt.10) then
        vfile(5:5)=char(id2+48)
      else
        vfile(5:6)=char(id1+48)//char(id2+48)
      end if
c
      open(13, file=vfile, form='unformatted', status='old')
      open(14, file='vel.mod', form='unformatted', status='old')
      open(35, file='for.header', status='old')
c
c     read in model dimensions
c
      read(35,115) xmin,xmax,ymin,ymax,zmin,zmax,size,nx,ny,nz
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
           read(13) (vel(i),i=1,nx)
10         write(14) (vel(i),i=1,nx)
c
      open(38, file='log.file',status='old',form='formatted')
155   read(38,55,end=98) alog
55    format(a1)
      go to 155
c
98    backspace(38)
      write(io,45) imin,chimin,lambmin
      write(38,45) imin,chimin,lambmin
45    format(/'iteration of best model: ',i3/
     +        'minimum chi and lambda: ',2e15.6)
c
      open(15, file='lambda', status='old')
      open(16, file='select.in', status='old')
c
      read(16,selpar)
c
      lambda=lambmin*strfac
c
      write(15,75) -lambda
75    format(e15.6)
c
      if(chione.gt.0) then
        open(48, file='chi.one', status='old')
        read(48,*) chione
        rewind(48)
        write(48,*) chione+1
        close(48)
      end if
c
      stop
c
      end
