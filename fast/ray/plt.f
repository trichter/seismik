c
c     version 1.0  Nov 1993
c
c     plot routines for RAY
c                 
c     ----------------------------------------------------------------
c                 
      subroutine pltsrcbox(source,iscol,ixy,ixz,i3d)
c                 
c     plot source box
c                 
      include 'ray.par'
      include 'ray.com'
c 
      integer source(6)
      real xpl(5),ypl(5)
c
      call pcolor(iscol)
c
      xpl(1)=size*float(source(1)-1)/xscale+xzxorig
      ypl(1)=zmm+size*float(source(5)-1)/zscale+xzyorig
      xpl(2)=size*float(source(2)-1)/xscale+xzxorig
      ypl(2)=ypl(1)
      xpl(3)=xpl(2)
      ypl(3)=zmm+size*float(source(6)-1)/zscale+xzyorig
      xpl(4)=xpl(1)
      ypl(4)=ypl(3)
      xpl(5)=xpl(1)
      ypl(5)=ypl(1)
      if(ixz.eq.1) call line(xpl,ypl,5)
      ypl(1)=ymm+size*float(source(3)-1)/yscale+xyyorig
      ypl(2)=ypl(1)
      ypl(3)=ymm+size*float(source(4)-1)/yscale+xyyorig
      ypl(4)=ypl(3)
      ypl(5)=ypl(1)
      if(ixy.eq.1) call line(xpl,ypl,5)
c
      call empty
c
      return
      end         
c                 
c     ----------------------------------------------------------------
c
      subroutine plotxy
c                 
      include 'ray.par'
      include 'ray.com'
c
      ipen=3
      call pcolor(3)
      open(37, file='ellipse')
      read(37,1) 
1     format(' ')
100   read(37,*,end=998) x,y
      x=(x-xmin)/xscale+xzxorig
      y=(180.-y-ymax)/yscale+xyyorig
      call plot(x,y,ipen)
      ipen=2
      go to 100
c
998   ipen=3
      call pcolor(3)
      open(38, file='rect')
      read(38,1)
200   read(38,*,end=999) x,z
      x=(x-xmin)/xscale+xzxorig
      z=(z-zmax)/zscale+xzyorig
      call plot(x,z,ipen)
      ipen=2
      go to 200
c
999   return
      end
