c
c     version 1.1  Apr 1995
c                 
c     Plotting routines for FD
c
c     ----------------------------------------------------------------
c                 
      subroutine pltmodxz(igrid)
c                 
c     plot the x-z slice of the velocity model 
c                 
      include 'fd.par'
      include 'fd.com'
c
      real xcp(2),zcp(2) 
c
      if(iplots.eq.0) then 
        call plots(xwndow,ywndow,iroute)
        call segmnt(1)
        iplots=1
      end if
      call erase
      call pcolor(ifcol)
      if(iaxlab.eq.1) then
        call axtick(xmin,xmax,xtmin,xtmax,ntickx,ndecix)
        call axis(xzxorig,xzyorig,xmin,xmax,xmm,xscale,0.,1,
     +       xtmin,xtmax,ntickx,ndecix,'X (km)',6,albht)
        call axtick(zmin,zmax,ztmin,ztmax,ntickz,ndeciz)
        call axis(xzxorig,xzyorig,zmin,zmax,zmm,zscale,90.,1,
     +       ztmin,ztmax,ntickz,ndeciz,'Z (km)',6,albht) 
      end if      
c                 
      call pcolor(icol(1))
c
      if(igrid.eq.1) then
        zcp(1)=zmm+xzyorig
        zcp(2)=xzyorig
        do 10 i=1,nx
           xcp(1)=float(i-1)*size/xscale+xzxorig
           xcp(2)=xcp(1)
           call line(xcp,zcp,2)
10      continue
c
        xcp(1)=xzxorig
        xcp(2)=xmm+xzxorig
        do 20 i=1,nz
           zcp(1)=(zmin+float(i-1)*size-zmax)/zscale+xzyorig
           zcp(2)=zcp(1)
           call line(xcp,zcp,2)
20      continue
      end if
c
      call pcolor(icol(3))
      call box(xzxorig,xzyorig,xzxorig+xmm,xzyorig+zmm)     
c
      call pcolor(ifcol)
c
      call empty
c
      return      
      end         
c
c     ----------------------------------------------------------------
c                 
      subroutine pltmodyz(imodxz,igrid)
c                 
c     plot the y-z slice of the velocity model 
c                 
      include 'fd.par'
      include 'fd.com'
c
      real ycp(2),zcp(2) 
c
      if(iplots.eq.0) then 
        call plots(xwndow,ywndow,iroute)
        call segmnt(1)
        iplots=1
        call erase
      end if
      call pcolor(ifcol)
      if(iaxlab.eq.1) then
        call axtick(ymin,ymax,ytmin,ytmax,nticky,ndeciy)
        call axis(yzxorig,yzyorig,ymin,ymax,ymm,yscale,0.,1,
     +       ytmin,ytmax,nticky,ndeciy,'Y (km)',6,albht)
        call axtick(zmin,zmax,ztmin,ztmax,ntickz,ndeciz)
        call axis(yzxorig,yzyorig,zmin,zmax,zmm,zscale,90.,1,
     +       ztmin,ztmax,ntickz,ndeciz,'Z (km)',6,albht)
      end if
c
      call pcolor(icol(1))
c
      if(igrid.eq.1) then
        zcp(1)=zmm+yzyorig
        zcp(2)=yzyorig
        do 10 i=1,ny
           ycp(1)=float(i-1)*size/yscale+yzxorig
           ycp(2)=ycp(1)
           call line(ycp,zcp,2)
10      continue
c
        ycp(1)=yzxorig
        ycp(2)=xmm+yzxorig
        do 20 i=1,nz
           zcp(1)=(zmin+float(i-1)*size-zmax)/zscale+yzyorig
           zcp(2)=zcp(1)
           call line(ycp,zcp,2)
20      continue
      end if
c
      call pcolor(icol(3))
      call box(yzxorig,yzyorig,yzxorig+ymm,yzyorig+zmm)
c
      call pcolor(ifcol)
c
      call empty
c
      return
      end
c
c     ----------------------------------------------------------------
c                 
      subroutine pltmodxy(imodxz,imodyz,igrid)
c                 
c     plot the x-y slice of the velocity model 
c                 
      include 'fd.par'
      include 'fd.com'
c
      real xcp(2),ycp(2) 
c
      if(iplots.eq.0) then 
        call plots(xwndow,ywndow,iroute)
        call segmnt(1)
        iplots=1
        call erase
      end if
      if(iaxlab.eq.1) then
        if(imodxz.eq.0) then
          call axtick(xmin,xmax,xtmin,xtmax,ntickx,ndecix)
          call axis(xzxorig,xzyorig,xmin,xmax,xmm,xscale,0.,1,
     +         xtmin,xtmax,ntickx,ndecix,'X (km)',6,albht)
        end if
        call axtick(ymin,ymax,ytmin,ytmax,nticky,ndeciy)
        call axis(xyxorig,xyyorig,ymin,ymax,ymm,-yscale,90.,1,
     +       ytmin,ytmax,nticky,ndeciy,'Y (km)',6,albht) 
      end if      
c                 
      call pcolor(icol(1))
c
      if(igrid.eq.1) then
        ycp(1)=ymm+xyyorig
        ycp(2)=xyyorig
        do 10 i=1,nx
           xcp(1)=float(i-1)*size/xscale+xyxorig
           xcp(2)=xcp(1)
           call line(xcp,ycp,2)
10      continue
c
        xcp(1)=xyxorig
        xcp(2)=xmm+xyxorig
        do 20 i=1,ny
           ycp(1)=-(ymin+float(i-1)*size-ymax)/yscale+xyyorig
           ycp(2)=ycp(1)
           call line(xcp,ycp,2)
20      continue
      end if
c
      call pcolor(icol(3))
      call box(xyxorig,xyyorig,xyxorig+xmm,xyyorig+ymm)     
c
      call pcolor(ifcol)
c
      call empty
c
      return      
      end         
c
c     ----------------------------------------------------------------
c                 
      subroutine pltsource(xs,ys,zs,imodxz,imodxy,imodyz)
c                 
c     plot the source location in the model 
c                 
      include 'fd.par'
      include 'fd.com'
c
      if(imodxz.eq.0.and.imodxy.eq.0.and.imodyz.eq.0) return
c
      call pcolor(icol(2))
c
      if(imodxz.eq.1) then
        xsp=(xs-xmin)/xscale+xzxorig
        zsp=(zs-zmax)/zscale+xzyorig
        call ssymbl(xsp,zsp,2.*nodeht,3)
      end if
c
      if(imodyz.eq.1) then
        ysp=(ys-ymin)/yscale+yzxorig
        zsp=(zs-zmax)/zscale+yzyorig
        call ssymbl(ysp,zsp,2.*nodeht,3)
      end if
c
      if(imodxy.eq.1) then
        xsp=(xs-xmin)/xscale+xyxorig
        ysp=-(ys-ymax)/yscale+xyyorig
        call ssymbl(xsp,ysp,2.*nodeht,3)
      end if
c
      call pcolor(ifcol)
c
      call empty
c
      return
      end
c
c     ----------------------------------------------------------------
c                 
      subroutine pltsrcnode(ixs,iys,izs,imodxz,imodxy,imodyz)
c                 
c     plot the model node nearest the source location
c                 
      include 'fd.par'
      include 'fd.com'
c
      if(imodxz.eq.0.and.imodxy.eq.0.and.imodyz.eq.0) return
c
      xs=float(ixs-1)*size+xmin
      ys=float(iys-1)*size+ymin
      zs=float(izs-1)*size+zmin
c
      call pcolor(icol(4))
c
      if(imodxz.eq.1) then
        xsp=(xs-xmin)/xscale+xzxorig
        zsp=(zs-zmax)/zscale+xzyorig
        call ssymbl(xsp,zsp,1.5*nodeht,4)
      end if
c
      if(imodyz.eq.1) then
        ysp=(ys-ymin)/yscale+yzxorig
        zsp=(zs-zmax)/zscale+yzyorig
        call ssymbl(ysp,zsp,1.5*nodeht,4)
      end if
c
      if(imodxy.eq.1) then
        xsp=(xs-xmin)/xscale+xyxorig
        ysp=-(ys-ymax)/yscale+xyyorig
        call ssymbl(xsp,ysp,1.5*nodeht,4)
      end if
c
      call pcolor(ifcol)
c
      call empty
c
      return
      end
c
c     ----------------------------------------------------------------
c                 
      subroutine plotnode(ix,iy,iz,imodxz,imodxy,imodyz)
c                 
c     plot a single model node
c                 
      include 'fd.par'
      include 'fd.com'
c
      xp=float(ix-1)*size+xmin
      yp=float(iy-1)*size+ymin
      zp=float(iz-1)*size+zmin
c
      if(imodxz.eq.1) then
        if(icflag.gt.0.or.xzplot(ix,iz).eq.0.or.
     +  (xzplot(ix,iz).ne.2.and.ipltnode.eq.8).or.
     +  (xzplot(ix,iz).ne.3.and.xzplot(ix,iz).ne.2.
     +  and.ipltnode.eq.6)) then
          xsp=(xp-xmin)/xscale+xzxorig
          zsp=(zp-zmax)/zscale+xzyorig
          call dot(xsp-n2,zsp+n2,nodeht,icol(ipltnode))
          if(ipltnode.eq.5.or.ipltnode.eq.7) then
            xzplot(ix,iz)=1
          else
            if(ipltnode.eq.8) then
              xzplot(ix,iz)=2
            else
              xzplot(ix,iz)=3
            end if
          end if
          nplot=nplot+1
        end if
      end if
c 
      if(imodyz.eq.1) then
        if(icflag.gt.0.or.yzplot(iy,iz).eq.0.or.
     +  (yzplot(iy,iz).ne.2.and.ipltnode.eq.8).or.
     +  (yzplot(iy,iz).ne.3.and.yzplot(iy,iz).ne.2.
     +  and.ipltnode.eq.6)) then
          ysp=(yp-ymin)/yscale+yzxorig
          zsp=(zp-zmax)/zscale+yzyorig
          call dot(ysp-n2,zsp+n2,nodeht,icol(ipltnode))
          if(ipltnode.eq.5.or.ipltnode.eq.7) then
            yzplot(iy,iz)=1
          else
            if(ipltnode.eq.8) then 
              yzplot(iy,iz)=2
            else
              yzplot(iy,iz)=3
            end if
          end if
          nplot=nplot+1
        end if
      end if
c
      if(imodxy.eq.1) then
        if(icflag.gt.0.or.xyplot(ix,iy).eq.0.or.
     +  (xyplot(ix,iy).ne.2.and.ipltnode.eq.8).or.
     +  (xyplot(ix,iy).ne.3.and.xyplot(ix,iy).ne.2.
     +  and.ipltnode.eq.6)) then
          xsp=(xp-xmin)/xscale+xyxorig
          ysp=-(yp-ymax)/yscale+xyyorig
          call dot(xsp-n2,ysp+n2,nodeht,icol(ipltnode))
          if(ipltnode.eq.5.or.ipltnode.eq.7) then
            xyplot(ix,iy)=1
          else
            if(ipltnode.eq.8) then
              xyplot(ix,iy)=2
            else
              xyplot(ix,iy)=3
            end if
          end if
          nplot=nplot+1
        end if
      end if
c
      call empty
c
      return
      end
