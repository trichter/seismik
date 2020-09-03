c
c     version 1.0  Jul 1998
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |                  *****  R E S E T  ****                      |
c     |                                                              |
c     |       Assign velocities to all nodes above an interface      |
c     |    according to a prior model and maintain all velocities    |
c     |                   within specified bounds                    |
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
c        10 -- input:  input 3D data grid
c
c        11 -- output:  output 3D data grid
c
c        12 -- input: bathymetry
c
c        35 -- input:  parameters describing size of 3D grid
c
c
c     ----------------------------------------------------------------
c
c
      include 'ray.par'
c
      parameter(npair=10)
      character*72 filei,fileo,fileb
      integer*2 vel(nxmax,nymax,nzmax),velb(nxmax,nymax,nzmax)
      real bath(nxmax,nymax),
     +     v0(npair)/npair*-1./,vmin(npair),vmax(npair)
c
      namelist /rspar/ ibound,v0,vmin,vmax,filei,fileo,fileb,idump,
     +                 interface
c
      interface=1
      ibound=0
      idump=0
c
      write(6,335)
335   format(
     +'RESET: assign prior model velocities and bounds')
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
      open(14, file='reset.in',status='old')
c
      read(14,rspar)
c
      open(10, file=filei, form='unformatted',status='old')
      open(13, file=fileb, form='unformatted',status='old')
      open(11, file=fileo, form='unformatted')
      if(interface.eq.1) 
     +open(12, file='bathymetry', form='unformatted',status='old')
      open(35, file='for.header', status='old')
      if(idump.eq.1) open(69, file='i.out')
c
c     read in model dimensions
c
      read(35,115) xmin,xmax,ymin,ymax,zmin,zmax,size,nx,ny,nz
115   format(7f10.3,3i10)
c
      if(nx.gt.nxmax.or.ny.gt.nymax) then
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
      if(ibound.eq.1) then
        nv=0
        do i=1,npair
           if(v0(i).lt.0.) go to 1998
           v0(i)=1000.*v0(i)
           vmin(i)=1000.*vmin(i)
           vmax(i)=1000.*vmax(i)
           nv=nv+1
        enddo
1998    if(nv.eq.0) then
          write(0,*) '***  must specify v0, vmin and vmax values  ***'
          stop
        end if
      end if
c
      zinc=(zmax-zmin)/float(nz-1)
      nfix=0
c
      if(interface.eq.1) then
        do 10 j=1,ny
10         read(12) (bath(i,j),i=1,nx)
      end if
c
      do 20 k=1,nz
         do 20 j=1,ny
            read(10) (vel(i,j,k),i=1,nx)
20          read(13) (velb(i,j,k),i=1,nx)
c
      do 30 i=1,nx
         do 30 j=1,ny 
            if(interface.eq.1) then
              ns1=nint((bath(i,j)/1000.-zmin)/size)+1
              do 31 k=1,ns1
31               vel(i,j,k)=velb(i,j,k)
            end if
c
            if(ibound.eq.1.and.ns1.lt.nz) then
              do 32 k=ns1+1,nz
c
               iflag=0
               do l=2,nv
                  if(float(velb(i,j,k)).le.v0(l).or.l.eq.nv) then
                    vpmin=(vmin(l)-vmin(l-1))/(v0(l)-v0(l-1))*
     +                    (velb(i,j,k)-v0(l-1))+vmin(l-1)
                    vpmax=(vmax(l)-vmax(l-1))/(v0(l)-v0(l-1))*
     +                    (velb(i,j,k)-v0(l-1))+vmax(l-1)
                    if(float(vel(i,j,k)).lt.vpmin) then
                      iflag=1
                      vf=vpmin
                      go to 1999
                    end if
                    if(float(vel(i,j,k)).gt.vpmax) then
                      iflag=1
                      vf=vpmax
                    end if
                    go to 1999
                  end if
               enddo
c
1999           if(iflag.eq.1) then
                 nfix=nfix+1
                 if(idump.eq.1) 
     +             write(69,95) i,j,k,velb(i,j,k),vel(i,j,k),vf
95               format(3i6,2i10,f10.3)
                 vel(i,j,k)=nint(vf)
               end if

32            continue
            end if
30    continue
c
      do 40 k=1,nz
         do 40 j=1,ny
40          write(11) (vel(i,j,k),i=1,nx)
c
      open(38, file='log.file',status='old',form='formatted')
155   read(38,55,end=98) alog
55    format(a1)
      go to 155
c
98    backspace(38)
      write(38,45) nfix
45    format('number of velocities set to bounds: ',i10)
c
      stop
      end
