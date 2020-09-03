c
c     version 1.0  Nov 1993
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |            *****  D I F F E R E N C E  ****                  |
c     |                                                              |
c     |          Take the difference of two 3D data grids            |
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
c        12 -- output: difference (file1 - file2)
c
c        35 -- input:  parameters describing size of 3D grid
c
c
c     ----------------------------------------------------------------
c
c
      include 'ray.par'
c
      character*72 file1,file2
      integer*2 veli1(nxmax),veli2(nxmax),imaxdiff
c
85    format(a72)
      write(*, fmt="(/'Enter file name')")
      read(5,85) file2
      write(*, fmt="(/'Enter background file name')")
      read(5,85) file1
c
      open(10, file=file1, form='unformatted', status='old')
      open(11, file=file2, form='unformatted', status='old')
      open(12, file='diff.out', form='unformatted')
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
      imaxdiff=0
      nk=0
c 
      do 10 k=1,nz
         do 11 j=1,ny
           read(10,end=999) (veli1(i),i=1,nx)
           read(11,end=999) (veli2(i),i=1,nx)
           do 21 i=1,nx 
              veli2(i)=nint(float(veli2(i)-veli1(i))/
     +                 veli1(i)*10000.)
              if(abs(veli2(i)).gt.imaxdiff) imaxdiff=abs(veli2(i))
21         continue
           write(12) (veli2(i),i=1,nx)
11       continue
         nk=nk+1
10    continue
c
      go to 200
c
999   if(nk.eq.1) then
        write(6,155) nk
155     format(/'>>> data file contains ',i3,' surface'/)
      else
        write(6,156) nk
156     format(/'>>> data file contains ',i3,' surfaces'/)
      end if
c
200   if(imaxdiff.eq.0) then
        write(6,3)
3       format(/'***  the two files are the same  ***'/)
      else
        write(6,4) imaxdiff
4       format(/'maximum difference between the two files: ',i8/)
      end if
c
      stop
      end
