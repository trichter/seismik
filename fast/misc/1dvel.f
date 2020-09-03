c
c     version 1.0  Aug 1994
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |                  *****  1 D V E L  ****                      |
c     |                                                              |
c     |        calculate the average 1D model from a 3D model        |
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
c        11 -- output:  output 3D model (that is the 1D average model)
c
c        12 -- input: num.cell   
c
c        13 -- output:  1D model ascii file
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
      integer*2 vel(nxmax,nymax,nzmax),
     +          numnode(nxmax,nymax,nzmax),vave
c
      write(*, fmt="(/'Enter input file name')")
      read(5,85) filei
85    format(a72)
      write(*, fmt="(/'Enter output file name (default is input)')")
      read(5,85) fileo
      if(fileo.eq.'') fileo=filei
      write(*, fmt="(/'Enter  1  to use num.cell')")
      read(5,*) inum
c
      open(10, file=filei, form='unformatted', status='old')
      open(11, file=fileo, form='unformatted')
      open(13, file='1dvel.out')
      open(14, file='1dsig.out')
      open(15, file='1drel.out')
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
      if(inum.eq.1) then
        open(12, file='num.cell', form='unformatted', status='old')
        do 10 k=1,nz
           do 10 j=1,ny
10            read(12) (numnode(i,j,k),i=1,nx)
      end if
c
      do 20 k=1,nz
         do 20 j=1,ny
20          read(10) (vel(i,j,k),i=1,nx)
c
      fnxny=float(nx*ny)
      zinc=(zmax-zmin)/float(nz-1)
c
      if(inum.ne.1) then
        do 30 k=1,nz
           depth=zmin+zinc*float(k-1)
           sum=0.
           sum2=0.
           do 40 j=1,ny
              do 40 i=1,nx 
                 sum=sum+vel(i,j,k) 
40               sum2=sum2+float(vel(i,j,k))**2 
           fvave=sum/fnxny
           vave=nint(sum/fnxny)
           sigma=sqrt(abs(-fnxny*fvave**2+sum2)/(fnxny-1))
           write(13,45) -depth,vave/1000.
           write(14,45) -depth,sigma/1000.
           write(15,45) -depth,100.*sigma/vave
45         format(2f10.3)
c
           do 50 j=1,ny
              do 60 i=1,nx
60               vel(i,j,k)=vave
              write(11) (vel(i,j,k),i=1,nx)
50         continue
30      continue
      else
       do 70 k=1,nz
           depth=zmin+zinc*float(k-1)
           num=0
           sum=0.
           sum2=0.
           do 80 j=1,ny
              do 80 i=1,nx
                 if(numnode(i,j,k).gt.0) then
                   num=num+1
                   sum=sum+vel(i,j,k)
                   sum2=sum2+float(vel(i,j,k))**2
                 end if
80         continue
           if(num.gt.1) then
             fvave=sum/num
             vave=nint(sum/num)
             sigma=sqrt(abs(-num*fvave**2+sum2)/(num-1))
           else
             sum=0.
             do 140 j=1,ny
                do 140 i=1,nx 
140                sum=sum+vel(i,j,k) 
             fvave=sum/fnxny
             vave=nint(sum/fnxny)
             sigma=0.
           end if
           write(13,45) -depth,vave/1000.
           write(14,45) -depth,sigma/1000.
           write(15,45) -depth,100.*sigma/vave
c
           do 90 j=1,ny
              do 110 i=1,nx
110              vel(i,j,k)=vave
              write(11) (vel(i,j,k),i=1,nx)
90         continue
70      continue 
      end if
c
      stop
      end
