c      
c     version 1.1  Apr 1995
c
c     common blocks for FD
c           
c        patched for g77 compilation
c
c        Drew Brenders (brenders@geoladm.geol.queensu.ca)
c
c     ----------------------------------------------------------------
c
      integer icol(8),colour(ncolour),stencil(8),sidelimit(6),
     +        hwside(6),tminnode(6,3),sidestop(6)
      integer*2 xyplot(nxmax,nymax),side(6),yzplot(nymax,nzmax),
     +           inter(nxmax,nymax)
c
c     use next 2 lines for 3D models
c
c      integer*2 xzplot(nxmax,nzmax)
c      real time(0:nxmax+1,0:nymax+1,0:nzmax+1)
c
c     use next 2 lines for 2D models
c
      integer*2 xzplot(nxmax,nymax)
      real time(0:nxmax+1,nymax,0:nzmax+1)
c
      real tminside(6),vel(nxmax,nymax,nzmax),nodeht,n2
c
      common /blk1/ xmin,xmax,xtmin,xtmax,xmm,ndecix,xscale,ntickx,
     +          ymin,ymax,ytmin,ytmax,ymm,ndeciy,yscale,nticky,
     +          zmin,zmax,ztmin,ztmax,zmm,ndeciz,zscale,ntickz,
     +          nodeht,albht,iplots,xzxorig,xzyorig,xyxorig,
     +          xyyorig,sep,icol,ncol,colour,
     +          iaxlab,yzxorig,yzyorig,n2,icflag
      common /blk2/ stencil,size,size2,nx,ny,nz,nodestotal,
     +              sidelimit,nside,i2d,ntup,nstop,sidestop,inter
      common /blk3/ n10,n10up,noplot,xzplot,xyplot,yzplot
      common /blk4/ nplot,ncalc,nnegsqrt,ncaus,nstencil,
     +              nsrccalc,iside,ipltnode,inode,iwrite,ireverse,
     +              iclear,istencil,io,interface
      common /blk5/ side,tminside,tminnode,hwside,tptmin,tminnew
      common /blk6/ time
      common /blk7/ vel
      common /cplot/ iplot,isep,iseg,nseg,xwndow,ywndow,ibcol,ifcol,sf
c
