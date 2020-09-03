c
c     version 1.0  Nov 1993
c
c     Block data for RAY
c
c     patched for g77 compilation 
c     Drew Brenders (brenders@geoladm.geol.queensu.ca)
c                 
c     ----------------------------------------------------------------
c                 
      block data  
c                 
c     assign default and initial values to common block parameters
c                 
      include 'ray.par'
      include 'ray.com'
c                 
      data 
     +  xmm,ndecix,ntickx/250.,-2,-1/,
     +  ymm,ndeciy,nticky/75.,-2,-1/,
     +  zmm,ndeciz,ntickz/75.,-2,-1/,
     +  tmm,ndecit,ntickt/75.,-2,-1/,
     +  xomm,ndecixo,ntickxo/250.,-2,-1/,
     +  xtmin,xtmax,ytmin,ytmax,ztmin,ztmax,tttmin,tttmax,
     +  xotmin,xotmax/10*-999999./,
     +  albht/2.5/,
     +  iplot,iplots,sep,iseg,nseg/1,0,7.5,0,0/,
     +  xwndow,ywndow/2*0./,colour/ncolour*-1/,
     +  sf,ibcol,ifcol/1.2,0,1/
c                 
      end         
