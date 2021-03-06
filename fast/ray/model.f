c
c     version 1.0  Nov 1993
c
c     model related routines for RAY
c                 
c     ----------------------------------------------------------------
c                 
      subroutine find(xp,yp,zp,ic,jc,kc)
c                 
c     find which model cell the point (xp,yp,zp) is in
c                 
      include 'ray.par'
      include 'ray.com'
c
      ic=int((xp-xmin)/size)+1
      jc=int((yp-ymin)/size)+1
      kc=int((zp-zmin)/size)+1
c
      return
      end         
c
c     ----------------------------------------------------------------
c
      subroutine findnode(xp,yp,zp,ixp,iyp,izp)
c
c     find closest model node to the point (xp,yp,zp)
c
      include 'ray.par'
      include 'ray.com'
c
      ixp=nint((xp-xmin)/size)+1
      iyp=nint((yp-ymin)/size)+1
      izp=nint((zp-zmin)/size)+1
c
      return
      end
c                 
c     ----------------------------------------------------------------
c                 
      subroutine grad(xp,yp,zp,xg,yg,zg)
c                 
c     calculate normal vector to local traveltime field for 
c     current ray point
c                 
      include 'ray.par'
      include 'ray.com'
c
      ixp=int((xp-xmin)/size)+1
      iyp=int((yp-ymin)/size)+1
      izp=int((zp-zmin)/size)+1
c
      x1=xp-sized2
      x2=xp+sized2
      ix1=int((x1-xmin)/size)+1
      ix2=int((x2-xmin)/size)+1
c

      xg=timeinterp(x1,yp,zp,ix1,iyp,izp)-
     +   timeinterp(x2,yp,zp,ix2,iyp,izp)
c
      if(ny.ne.1) then
        y1=yp-sized2
        y2=yp+sized2
        iy1=int((y1-ymin)/size)+1
        iy2=int((y2-ymin)/size)+1
c

        yg=timeinterp(xp,y1,zp,ixp,iy1,izp)-
     +     timeinterp(xp,y2,zp,ixp,iy2,izp)
      else
        yg=0.
      end if
c
      z1=zp-sized2
      z2=zp+sized2
      iz1=int((z1-zmin)/size)+1
      iz2=int((z2-zmin)/size)+1
c

      zg=timeinterp(xp,yp,z1,ixp,iyp,iz1)-
     +   timeinterp(xp,yp,z2,ixp,iyp,iz2)
c
      return
      end         
c                 
c     ----------------------------------------------------------------
c                 
      subroutine intersect(ic,jc,kc,xs,ys,zs,ixs,iys,izs,xg,yg,zg,
     +              xp,yp,zp,xint,yint,zint,imode)
c                 
c     calculate intersection point of ray with cell boundary
c                 
      include 'ray.par'
      include 'ray.com'
c
      real :: rparameter(6)
      integer i,xxx,status_read
      logical file_exist
c
c      print*,'INTERSECT-1:',ic,jc,kc,xs,ys,zs,ixs,iys,izs,xg,yg,zg
c      print*,'INTERSECT-2:',xp,yp,zp,xint,yint,zint,imode

      xxx=0
      inquire(file='read_parameter.txt',exist=file_exist)
      if(file_exist)then
      open(unit=100000,file='read_parameter.txt',status='old',
     +form='formatted')
      read(100000,*,iostat=status_read)xxx
      if(status_read/=0)xxx=0
      close(100000)
      endif
1000  continue
c
      do 70 i=1,6
70       rparameter(i)=1.e21
c
      if(xg.lt.0.) then
c
c       try left side of cell
c
        if(ic.eq.ixs.and.imode.eq.0) then
          x1i=xs
        else
          x1i=xmodel(ic)
        end if
        rparameter(1)=(x1i-xp)/xg
      end if 
c 
      if(xg.gt.0.) then
c
c       try right side of cell
c
        if(ic.eq.ixs.and.imode.eq.0) then
          x2i=xs
        else
          x2i=xmodel(ic+1)
        end if
        rparameter(2)=(x2i-xp)/xg
      end if
c
      if(yg.lt.0.) then
c
c       try back side of cell
c
        if(jc.eq.iys.and.imode.eq.0) then
          y3i=ys
        else
          y3i=ymodel(jc)
        end if
        rparameter(3)=(y3i-yp)/yg
      end if 
c
      if(yg.gt.0.) then
c
c       try front side of cell
c
        if(jc.eq.iys.and.imode.eq.0) then
          y4i=ys
        else
          y4i=ymodel(jc+1)
        end if
        rparameter(4)=(y4i-yp)/yg
      end if
c
      if(zg.lt.0.) then
c
c       try top side of cell
c
        if(kc.eq.izs.and.imode.eq.0) then
          z5i=zs
        else
          z5i=zmodel(kc)
        end if
        rparameter(5)=(z5i-zp)/zg
      end if 
c
      if(zg.gt.0.) then
c
c       try bottom side of cell
c
        if(kc.eq.izs.and.imode.eq.0) then
          z6i=zs
        else
          z6i=zmodel(kc+1)
        end if
        rparameter(6)=(z6i-zp)/zg
      end if
c if rparamater is nearly zero (e.g. -8*10^-6) the following conditions are fail.
c thats the reason why we have to set rparameter to zero
      if(xxx==1)print*,'PARAMETER1:',rparameter
      do i=1,6
        if(sqrt(rparameter(i)**2)<0.00001)rparameter(i)=0.00001
      end do
      if(xxx==1)print*,'PARAMETER2:',rparameter
      pmin=minval(rparameter(1:6))
      if(xxx==1)print*,pmin
      ni=0
      xint=0.
      yint=0.
      zint=0.
c
      if(.999999*rparameter(1).le.pmin) then
c        print*,'IF1'
        xint=x1i
        yint=yg*rparameter(1)+yp
        zint=zg*rparameter(1)+zp
        ni=ni+1
        ic=ic-1
      else   
        if(.999999*rparameter(2).le.pmin) then
c        print*,'IF2'
          xint=x2i
          yint=yg*rparameter(2)+yp
          zint=zg*rparameter(2)+zp
          ni=ni+1
          ic=ic+1
        end if
      end if
      if(.999999*rparameter(3).le.pmin) then
c        print*,'IF3'
        xint=xint+xg*rparameter(3)+xp
        yint=yint+y3i
        zint=zint+zg*rparameter(3)+zp
        ni=ni+1
        jc=jc-1
      else   
        if(.999999*rparameter(4).le.pmin) then
c        print*,'IF4'
          xint=xint+xg*rparameter(4)+xp
          yint=yint+y4i
          zint=zint+zg*rparameter(4)+zp
          ni=ni+1
          jc=jc+1
        end if
      end if
      if(.999999*rparameter(5).le.pmin) then
c        print*,'IF5'
        xint=xint+xg*rparameter(5)+xp
        yint=yint+yg*rparameter(5)+yp
        zint=zint+z5i
        ni=ni+1
        kc=kc-1
      else   
        if(.999999*rparameter(6).le.pmin) then
c        print*,'IF6'
          xint=xint+xg*rparameter(6)+xp
          yint=yint+yg*rparameter(6)+yp
          zint=zint+z6i
          ni=ni+1
          kc=kc+1
        end if
      end if
c
      xint=xint/ni
      if(ny.ne.1) then
        yint=yint/ni
      else
        yint=yp
      end if
      zint=zint/ni
c      print*,xint,yint,zint,ni
c
      return
      end
c                 
c     ----------------------------------------------------------------
c                 
      subroutine straight(ic,jc,kc,xs,ys,zs,ixs,iys,izs,iflag)
c                 
c     calculate straight-line ray path inside source box to the source
c                 
      include 'ray.par'
      include 'ray.com'
c
      iflag=0
c
      xp=xray(npts) 
      yp=yray(npts) 
      zp=zray(npts) 
c
      xg=xs-xp
      if(ny.ne.1) then
        yg=ys-yp
      else
        yg=0.
      end if
      zg=zs-zp
c
1000  continue
c
      call intersect(ic,jc,kc,xs,ys,zs,ixs,iys,izs,xg,yg,zg,
     +               xp,yp,zp,xint,yint,zint,0)
c
      xp=xint
      yp=yint
      zp=zint
      npts=npts+1
      xray(npts)=xp
      yray(npts)=yp
      zray(npts)=zp
c
      cell(npts,1)=ic
      cell(npts,2)=jc
      cell(npts,3)=kc
c
      if(abs(xp-xs).lt.smin.and.abs(yp-ys).lt.smin.
     +   and.abs(zp-zs).lt.smin) go to 999
c
c
      if(npts.eq.nray) then
c
c       ray consists of too many points
c
        iflag=1
        go to 999
      end if
c
      if(ic.ne.ixs.or.jc.ne.iys.or.kc.ne.izs) go to 1000
c
      npts=npts+1
      xray(npts)=xs
      yray(npts)=ys
      zray(npts)=zs
c
      cell(npts,1)=ic
      cell(npts,2)=jc
      cell(npts,3)=kc
c       
999   return
      end
c
c     ----------------------------------------------------------------
c
      function bndinterp(xp,yp,ixc,iyc,itype)
c
c     calculate the boundary depth at the point (xp,yp,zp), in general not
c     a grid point
c
      include 'ray.par'
      include 'ray.com'
c
      if(itype.eq.1) then
        z1=boundary(ixc,iyc)
        z2=boundary(ixc+1,iyc)
        z3=boundary(ixc,iyc+1)
        z4=boundary(ixc+1,iyc+1)
      else
        z1=inter(ixc,iyc)
        z2=inter(ixc+1,iyc)
        z3=inter(ixc,iyc+1)
        z4=inter(ixc+1,iyc+1)
      end if
c
      xfrac=xp-(float(ixc-1)*size+xmin)
      yfrac=yp-(float(iyc-1)*size+ymin)
c
      zleft= (z3-z1)*yfrac/size+z1
      zright=(z4-z2)*yfrac/size+z2
c
      bndinterp=((size-xfrac)*zleft+xfrac*zright)/size
c
      return
      end
