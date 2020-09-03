c
c     version 1.0  Jul 1994
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |               ********  R E G R I D 2  *****                 |
c     |                                                              |
c     |   Resample 3D inverse cell mode onto 3D forward node model   |
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
c        10 -- input:  inverse 3D data grid
c
c        11 -- output:  forward 3D data grid
c
c        35 -- input:  parameters describing size of forward 3D grid
c
c        36 -- input:  parameters describing size of inverse 3D grid
c
c
c     ----------------------------------------------------------------
c
c
      include 'ray.par'
      common /blk/ xmin,ymin,zmin,xinc,yinc,zinc,nxi,nyi,nzi,veln
c
      real veli(nximax+1,nyimax+1,nzimax+1),velr(nxmax),
     +     veln(nximax+1,nyimax+1,nzimax+1)
      integer*2 ivel(nxmax)
      character*72 file1,file2
      character*1 r4
c
      write(6,335)
335   format('REGRID2: map inverse model to forward model')
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
      write(io, fmt="(/'Enter input file name')")
      read(5,85) file1
85    format(a72)
      write(io, fmt="(/'Enter output file name')")
      read(5,85) file2
      if(file2.eq.'') file2=file1
      write(io, fmt="(/'Enter  r  for r*4 file (default is i*2)')")
      read(5,95) r4
95    format(a1)
      write(io, fmt="(/'Enter number of smoothing applications')")
      read(5,*) nsmooth
      write(io, fmt="(/
     +  'Enter x,y,z smoothing operator half width')")
      read(5,*) nxsmooth,nysmooth,nzsmooth
c
      open(10, file=file1, form='unformatted', status='old')
      open(35, file='for.header', status='old')
      open(36, file='inv.header', status='old')
c
c     read in model dimensions
c
      read(35,115) xmin,xmax,ymin,ymax,zmin,zmax,size,nx,ny,nz
115   format(7f10.3,3i10)
c
      if(nx.gt.nxmax.or.ny.gt.nymax.or.nz.gt.nzmax) then
        write(0,9)
9       format(/'***  model is too large  ***'/)
        stop
      end if
c 
      nodestotal=nx*ny*nz 
c
      read(36,*) nxi,nyi,nzi
c
      if(nxi.lt.1.or.nyi.lt.1.or.nzi.lt.1.or.
     +   nxi.ge.nxmax.or.nyi.ge.nymax.or.nzi.ge.nzmax) then
        write(0,8)
8       format(/'***  original file has invalid size  ***'/)
        stop
      end if
c
      nodestotali=nxi*nyi*nzi
c
      write(io,25) xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,nodestotal 
25    format(/'model dimensions:'/ 
     +        '-----------------'/ 
     +        'xmin, xmax: ',2f10.3/'ymin, ymax: ',2f10.3/, 
     +        'zmin, zmax: ',2f10.3/ 
     +        'number of model nodes in each direction (x,y,z):',3i4/ 
     +        'total number of model nodes:                    ',i12) 
      write(io,35) size 
35    format('node spacing: ',f7.3,' km') 
      write(io,45) nxi,nyi,nzi,nodestotali
45    format(
     +'original number of cells in each direction (x,y,z):  ',3i4/
     +'original total number of cells:                      ',i12)
c
      xinc=(xmax-xmin)/float(nxi)
      yinc=(ymax-ymin)/float(nyi)
      zinc=(zmax-zmin)/float(nzi)
c
      if(r4.eq.'r') then
        do 10 k=1,nzi
           do 10 j=1,nyi
10            read(10) (veli(i,j,k),i=1,nxi)
      else
        do 20 k=1,nzi
           do 20 j=1,nyi
              read(10) (ivel(i),i=1,nxi)
              do 60 i=1,nxi
60               veli(i,j,k)=ivel(i)
20      continue
      end if
c
      close(10)
      open(11, file=file2, form='unformatted')
c
      do 30 k=1,nzi+1
         do 30 j=1,nyi+1
            do 300 i=1,nxi+1
               if(i.eq.1) then
                 imin=i
               else
                 imin=i-1
               end if
               if(j.eq.1) then
                 jmin=j
               else
                 jmin=j-1
               end if
               if(k.eq.1) then
                 kmin=k
               else
                 kmin=k-1
               end if
               if(i.eq.nxi+1) then
                 imax=i-1
               else
                 imax=i
               end if
               if(j.eq.nyi+1) then
                 jmax=j-1
               else
                 jmax=j
               end if
               if(k.eq.nzi+1) then
                 kmax=k-1
               else
                 kmax=k
               end if
               veln(i,j,k)=0.
               nc=0
               do 400 ii=imin,imax
                  do 400 jj=jmin,jmax
                     do 400 kk=kmin,kmax
                        nc=nc+1
                        veln(i,j,k)=veln(i,j,k)+veli(ii,jj,kk)
400            continue
               veln(i,j,k)=veln(i,j,k)/float(nc)
300         continue
30    continue
c
      if(nsmooth.gt.0) then
      do 180 l=1,nsmooth
         do 1110 k=1,nzi+1
            do 1110 j=1,nyi+1
               do 1110 i=1,nxi+1
                  imin=max(1,i-nxsmooth)
                  imax=min(nxi+1,i+nxsmooth)
                  jmin=max(1,j-nysmooth)
                  jmax=min(nyi+1,j+nysmooth)
                  kmin=max(1,k-nzsmooth)
                  kmax=min(nzi+1,k+nzsmooth)
                  if(i.eq.1) then
                    sum=0.
                    do 1120 ii=imin,imax
                       do 1120 jj=jmin,jmax
                          do 1120 kk=kmin,kmax
1120                         sum=sum+veln(ii,jj,kk)
                  else
                    if(imin.gt.iminl) then
                      do 1130 jj=jmin,jmax
                         do 1130 kk=kmin,kmax
1130                        sum=sum-veln(imin-1,jj,kk)
                    end if
                    if(imax.gt.imaxl) then
                      do 1140 jj=jmin,jmax
                         do 1140 kk=kmin,kmax
1140                        sum=sum+veln(imax,jj,kk)
                    end if
                  end if
                  iminl=imin
                  imaxl=imax
                  nsum=(imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1)
1110              veli(i,j,k)=sum/float(nsum)
c
         do 170 i=1,nxi+1
            do 170 j=1,nyi+1
               do 170 k=1,nzi+1
170               veln(i,j,k)=veli(i,j,k)
c
         write(io,75) l
75       format('smoothing iteration: ',i5)
c
180   continue
      end if
c
      if(r4.eq.'r') then
        do 40 k=1,nz
           zp=float(k-1)*size+zmin
           do 40 j=1,ny
              yp=float(j-1)*size+ymin
              do 50 i=1,nx
                 xp=float(i-1)*size+xmin
                 velr(i)=vinterp(xp,yp,zp)
50            continue
              write(11) (velr(i),i=1,nx)
40      continue
      else
        do 70 k=1,nz
           zp=float(k-1)*size+zmin
           do 70 j=1,ny
              yp=float(j-1)*size+ymin
              do 80 i=1,nx
                 xp=float(i-1)*size+xmin
                 ivel(i)=vinterp(xp,yp,zp)
80            continue
              write(11) (ivel(i),i=1,nx)
70      continue
      end if
c
      stop
      end
c
c     ----------------------------------------------------------------
c
      function vinterp(xp,yp,zp)
c
c     calculate the value at the point (xp,yp,zp), in general not
c     a grid point
c
      include 'ray.par'
      common /blk/ xmin,ymin,zmin,xinc,yinc,zinc,nxi,nyi,nzi,veln
      real veln(nximax+1,nyimax+1,nzimax+1)
c
      ixc=min(nxi,int((xp-xmin)/xinc)+1)
      iyc=min(nyi,int((yp-ymin)/yinc)+1)
      izc=min(nzi,int((zp-zmin)/zinc)+1)
c
      v1=veln(ixc,iyc,izc)
      v2=veln(ixc+1,iyc,izc)
      v3=veln(ixc,iyc+1,izc)
      v4=veln(ixc+1,iyc+1,izc)
      v5=veln(ixc,iyc,izc+1)
      v6=veln(ixc+1,iyc,izc+1)
      v7=veln(ixc,iyc+1,izc+1)
      v8=veln(ixc+1,iyc+1,izc+1)
      xfrac=xp-(float(ixc-1)*xinc+xmin)
      yfrac=yp-(float(iyc-1)*yinc+ymin)
      zfrac=zp-(float(izc-1)*zinc+zmin)
c
      vleft= (zinc*(v3-v1)*yfrac+yinc*(v5-v1)*zfrac+
     +       (v1-v3+v7-v5)*yfrac*zfrac)/(yinc*zinc)+v1
      vright=(zinc*(v4-v2)*yfrac+yinc*(v6-v2)*zfrac+
     +       (v2-v4+v8-v6)*yfrac*zfrac)/(yinc*zinc)+v2
c
      vinterp=((xinc-xfrac)*vleft+xfrac*vright)/xinc
c
      return
      end
