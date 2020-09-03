c
c     version 1.1  Apr 1995
c
c     Block data for FD
c                 
c     ----------------------------------------------------------------
c                 
      block data  
c                 
c     assign default and initial values to common block parameters
c                 
      include 'fd.par'
      include 'fd.com'
c                 
      data 
     +  xmin,xmax,xmm,ndecix,ntickx/0.,-99999.,250.,-2,-1/,
     +  ymin,ymax,ymm,ndeciy,nticky/0.,50.,75.,-2,-1/,
     +  zmin,zmax,zmm,ndeciz,ntickz/0.,50.,75.,-2,-1/,
     +  xtmin,xtmax,ytmin,ytmax,ztmin,ztmax/
     +  6*-999999./,nodeht,albht/.5,2.5/,
     +  iplot,iplots,sep,iseg,nseg/1,0,7.5,0,0/,
     +  xwndow,ywndow/2*0./,colour/ncolour*-1/,icol/8*-1/,
     +  sf,ibcol,ifcol/1.2,0,1/
c                 
      end         
