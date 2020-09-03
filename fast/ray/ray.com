c      
c     version 1.0  Jan 1994
c
c     common blocks for RAY
c 
c      patched for g77 compilation
c      
c      Drew Brenders (brenders@geoladm.geol.queensu.ca)
c           
c     ----------------------------------------------------------------
c
      integer colour(ncolour),cell(nray,3),row(nwrite),column(nwrite)
      integer*2 time(nxmax,nymax,nzmax),numray(nximax,nyimax,nzimax),
     +          one
      real xray(nray),yray(nray),zray(nray),
     +     xmodel(nxmax),ymodel(nymax),zmodel(nzmax),
     +     xplot(nray),yplot(nray),zplot(nray),lseg(nray),
     +   pert1(nximax,nyimax,nzimax),pert2(nximax,nyimax,nzimax),
     +   boundary(nxmax,nymax),dkernel(nwrite),
     +   imodel(nximax,nyimax,nzimax),inter(nxmax,nymax)
c
      common /blk1/ xmin,xmax,xtmin,xtmax,xmm,ndecix,xscale,ntickx,
     +     ymin,ymax,ytmin,ytmax,ymm,ndeciy,yscale,nticky,
     +     zmin,zmax,ztmin,ztmax,zmm,ndeciz,zscale,ntickz,
     +     ttmin,ttmax,tttmin,tttmax,tmm,ndecit,tscale,ntickt,
     +     xomin,xomax,xotmin,xotmax,xomm,ndecixo,xoscale,ntickxo,
     +     albht,iplots,xzxorig,xzyorig,xyxorig,xyyorig,xtxorig,
     +     xtyorig,sep,ncol,colour,iaxlab
c
      common /blk2/ size,size2,sizex2,sizex4,nx,ny,nz,boundary,
     +       xmodel,ymodel,zmodel,txsize,txsizex2,sized2,inx,iny,inz,
     +       nm,xii,yii,zii,inter
c
      common /blk3/ xray,yray,zray,xplot,yplot,zplot
c
      common /blk4/ lseg,nlc,istype,
     +              cell,ncount,npts,smin,imethod,nk,nl,one
c
      common /blk5/ time 
c
      common /blk6/ numray
c
      common /blk7/ pert1
c
      common /blk8/ pert2
c
      common /blk9/ imodel
c
      common /cplot/ iplot,isep,iseg,nseg,xwndow,ywndow,ibcol,ifcol,sf
c
