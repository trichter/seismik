c      
c     version 1.0  Nov 1993
c
c     common blocks for ZSLICE
c           
c     ----------------------------------------------------------------
c
      integer colour(ncolour)
      real*4 time(nxmax,nymax,nzmax),tplot(nx2dmax,ny2dmax)
c
      common /blk1/ xmin,xmax,xtmin,xtmax,xmm,ndecix,xscale,ntickx,
     +          ymin,ymax,ytmin,ytmax,ymm,ndeciy,yscale,nticky,
     +          zmin,zmax,ztmin,ztmax,zmm,ndeciz,zscale,ntickz,
     +          albht,iplots,xorig,yorig,ncol,colour,iaxlab
      common /blk2/ time,tplot,size,size2,nx,ny,nz,undef
      common /cplot/ iplot,isep,iseg,nseg,xwndow,ywndow,ibcol,ifcol,sf
c
