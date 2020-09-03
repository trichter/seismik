c
c     version 1.1  Apr 1995
c
c     time routine for FD
c
c     ----------------------------------------------------------------
c
      subroutine timenearsrc(xs,ys,zs,ixs,iys,izs,nptsrc,
     +        imodxz,imodxy,imodyz,istop,inear,vabove,vbelow)
c
c     time all points inside an nptsrc-sided box centred 
c     around the node nearest the source
c
      include 'fd.par'
      include 'fd.com'
c
c     first calculate velocity at the source
c
      call find(xs,ys,zs,ixc,iyc,izc)
      vs=velocity(xs,ys,zs,ixc,iyc,izc)
      if(iwrite.eq.1) write(io,1) vs
1     format('          velocity at source: ',f7.2)
c
c     now calculate time to all points within box around source
c
      il=(nptsrc-1)/2
c
      if(ixs-il.lt.1) then
        ixmin=1
      else
        ixmin=ixs-il
      end if
      side(1)=ixmin
      if(ixs+il.gt.nx) then
        ixmax=nx
      else
        ixmax=ixs+il
      end if
      side(2)=ixmax
      if(iys-il.lt.1) then
        iymin=1
      else
        iymin=iys-il
      end if
      side(3)=iymin
      if(iys+il.gt.ny) then
        iymax=ny
      else
        iymax=iys+il
      end if
      side(4)=iymax
      if(izs-il.lt.1) then
        izmin=1
      else
        izmin=izs-il
      end if
      side(5)=izmin
      if(izs+il.gt.nz) then
        izmax=nz
      else
        izmax=izs+il
      end if
      side(6)=izmax
c
      if(i2d.eq.1) then
        side(3)=iys
        side(4)=iys
      end if
c
      ipltnode=7
      do 10 i=side(1),side(2)
         do 10 j=side(3),side(4)
            do 10 k=side(5),side(6)
               if(inear.ne.1) then
                 vss=vs
               else
                 zpoint=float(k-1)*size+zmin
                 if(zpoint.lt.zs) then
                   vss=vabove
                 else
                   vss=vbelow
                 end if
               end if
               time(i,j,k)=timelinear(0.,xs,ys,zs,vss,i,j,k)
               ncalc=ncalc+1
               if(noplot.eq.0) call plotnode(i,j,k,imodxz,imodxy,imodyz)
               if(i.eq.side(1).and.time(i,j,k).lt.tminside(1)) then 
                 tminside(1)=time(i,j,k)
                 tminnode(1,1)=i
                 tminnode(1,2)=j
                 tminnode(1,3)=k
               end if
               if(i.eq.side(2).and.time(i,j,k).lt.tminside(2)) then
                 tminside(2)=time(i,j,k)
                 tminnode(2,1)=i
                 tminnode(2,2)=j
                 tminnode(2,3)=k
               end if
               if(j.eq.side(3).and.time(i,j,k).lt.tminside(3)) then
                 tminside(3)=time(i,j,k)
                 tminnode(3,1)=i
                 tminnode(3,2)=j
                 tminnode(3,3)=k
               end if
               if(j.eq.side(4).and.time(i,j,k).lt.tminside(4)) then
                 tminside(4)=time(i,j,k)
                 tminnode(4,1)=i
                 tminnode(4,2)=j
                 tminnode(4,3)=k
               end if
               if(k.eq.side(5).and.time(i,j,k).lt.tminside(5)) then
                 tminside(5)=time(i,j,k)
                 tminnode(5,1)=i
                 tminnode(5,2)=j
                 tminnode(5,3)=k
               end if
               if(k.eq.side(6).and.time(i,j,k).lt.tminside(6)) then
                 tminside(6)=time(i,j,k)
                 tminnode(6,1)=i
                 tminnode(6,2)=j
                 tminnode(6,3)=k
               end if
10    continue
c
      do 20 i=1,6
         if(side(i).eq.sidelimit(i)) nside=nside+1
20    continue
      if(istop.eq.1) then
        if(side(1).le.sidestop(1)) nstop=nstop+1
        if(side(2).ge.sidestop(2)) nstop=nstop+1
        if(side(3).le.sidestop(3)) nstop=nstop+1
        if(side(4).ge.sidestop(4)) nstop=nstop+1
        if(side(5).le.sidestop(5)) nstop=nstop+1
        if(side(6).ge.sidestop(6)) nstop=nstop+1
      end if
      nsrccalc=ncalc
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine readsrc(is,imodxz,imodxy,imodyz,istop,tmult)
c
c     read in the times of a pre-existing 3D time field
c
      include 'fd.par'
      include 'fd.com'
c
      character*11 sfile
      integer*2 itime(nxmax)
c
      data sfile/'fd  .source'/
c
      id1=is/10
      id2=is-id1*10
      sfile(3:4)=char(id1+48)//char(id2+48)
      open(29, file=sfile, form='unformatted', status='old')
      open(49, file='int.kmax',status='old')
c
      read(49,*) kmax
      close(49)
c
      if(iwrite.eq.1) write(io,1)
1     format('          source times input from file')
c
      side(1)=1
      side(2)=nx
      side(3)=1
      side(4)=ny
      side(5)=kmax
      side(6)=nz
c
c     read times to all points
c
      do 10 k=1,nz
         do 10 j=1,ny 
            read(29) (itime(i),i=1,nx)
            do 20 i=1,nx
               time(i,j,k)=itime(i)/tmult
               if(itime(i).eq.0) time(i,j,k)=1.e10
20          continue
10    continue
c
      nside=5
      nsrccalc=0
c
      return
      end
c
c     ----------------------------------------------------------------
c
      function timelinear(t1,x1,y1,z1,v1,ix2,iy2,iz2)
c
c     calculate the time from point 1 to point 2 assuming a straight
c     ray path and linear velocity gradient from velocity at point 1
c     to velocity at point 2
c
      include 'fd.par'
      include 'fd.com'
c
      x2=float(ix2-1)*size+xmin
      y2=float(iy2-1)*size+ymin
      z2=float(iz2-1)*size+zmin
c
      v2=size*vel(ix2,iy2,iz2)
      dist=((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)**.5
c
      if(abs(v1-v2).lt..01) then
        timelinear=t1+2.*dist/(v1+v2)
      else
        timelinear=t1+dist*log(v2/v1)/(v2-v1)
      end if
c
      return
      end
c
c     ----------------------------------------------------------------
c
      function timeinterp(xp,yp,zp,ixc,iyc,izc)
c
c     calculate the traveltime at the point (xp,yp,zp), in general not
c     a grid point
c
      include 'fd.par'
      include 'fd.com'
c
      t1=time(ixc,iyc,izc)
      t2=time(ixc+1,iyc,izc)
      t3=time(ixc,iyc+1,izc)
      t4=time(ixc+1,iyc+1,izc)
      t5=time(ixc,iyc,izc+1)
      t6=time(ixc+1,iyc,izc+1)
      t7=time(ixc,iyc+1,izc+1)
      t8=time(ixc+1,iyc+1,izc+1)
      xfrac=xp-(float(ixc-1)*size+xmin)
      yfrac=yp-(float(iyc-1)*size+ymin)
      zfrac=zp-(float(izc-1)*size+zmin)
c
      tleft= (size*(t3-t1)*yfrac+size*(t5-t1)*zfrac+
     +       (t1-t3+t7-t5)*yfrac*zfrac)/size2+t1 
      tright=(size*(t4-t2)*yfrac+size*(t6-t2)*zfrac+
     +       (t2-t4+t8-t6)*yfrac*zfrac)/size2+t2 
c
      timeinterp=((size-xfrac)*tleft+xfrac*tright)/size
c
      return
      end
c
c     ----------------------------------------------------------------
c
      function timeinterp2d(xp,zp,ixc,izc,iys)
c
c     calculate the traveltime at the point (xp,yp,zp), in general not
c     a grid point, for 2D time field
c
      include 'fd.par'
      include 'fd.com'
c
      t1=time(ixc,iys,izc)
      t2=time(ixc+1,iys,izc)
      t5=time(ixc,iys,izc+1)
      t6=time(ixc+1,iys,izc+1)
      xfrac=xp-(float(ixc-1)*size+xmin)
      zfrac=zp-(float(izc-1)*size+zmin)
c
      tleft= (t5-t1)*zfrac/size+t1 
      tright=(t6-t2)*zfrac/size+t2 
c
      timeinterp2d=((size-xfrac)*tleft+xfrac*tright)/size
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine timea(ix0,iy0,iz0,ix7,iy7,iz7)
c
c     calculate the time at point 7 using stencil A centred at point 0
c
      include 'fd.par'
      include 'fd.com'
c
      save=1./vaverage(ix0,iy0,iz0,ix7,iy7,iz7)
c
      t0=time(ix0,iy0,iz0)
      t1=time(ix7,iy0,iz0)
      t2=time(ix0,iy7,iz0)
      t3=time(ix7,iy7,iz0)
      t4=time(ix0,iy0,iz7)
      t5=time(ix7,iy0,iz7)
      t6=time(ix0,iy7,iz7)
c
      t12=t1-t2
      t24=t2-t4
      t41=t4-t1
      t35=t3-t5
      t56=t5-t6
      t63=t6-t3 
c
      sqrtterm=6.*save*save-(t12*t12+t24*t24+t41*t41+
     +                       t35*t35+t56*t56+t63*t63)
c
      if(sqrtterm.lt.0.) then          
        nnegsqrt=nnegsqrt+1
      else
        ta=t0+sqrt(.5*sqrtterm)
        if(ta.lt.max(t1,t2,t3,t4,t5,t6)) then
          ncaus=ncaus+1
        else
          if(ta.lt.tptmin) then
            tptmin=ta
            istencil=1
          end if
        end if
      end if
c
      nstencil=nstencil+1
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine timeb(ix0,iy0,iz0,ix1,iy1,iz1,ix2,iy2,iz2,
     +               ix3,iy3,iz3,ix4,iy4,iz4,ix5,iy5,iz5)
c
c     calculate the time at point 5 using stencil B and points 
c     0 through 4
c
      include 'fd.par'
      include 'fd.com'
c
      save=2./(vaverage(ix0,iy0,iz0,ix5,iy5,iz5)+
     +         vaverage(ix3,iy3,iz3,ix5,iy5,iz5))
c
      t0=time(ix0,iy0,iz0)
      t1=time(ix1,iy1,iz1)
      t2=time(ix2,iy2,iz2)
      t3=time(ix3,iy3,iz3)
      t4=time(ix4,iy4,iz4)
c
      t03=t0-t3
      t24=t2-t4
c
      sqrtterm=2.*save*save-(0.5*t03*t03+t24*t24)
c
      if(sqrtterm.lt.0.) then          
        nnegsqrt=nnegsqrt+1
      else
        tb=t1+sqrt(sqrtterm)
        if(tb.lt.max(t0,t2,t3,t4)) then
          ncaus=ncaus+1 
        else
          if(tb.lt.tptmin) then   
            tptmin=tb
            istencil=2
          end if
        end if
      end if
c
      nstencil=nstencil+1
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine timec(ix0,iy0,iz0,ix1,iy1,iz1,ix2,iy2,iz2,
     +               ix3,iy3,iz3,ix4,iy4,iz4,ix5,iy5,iz5,
     +               ix6,iy6,iz6,ix7,iy7,iz7,
     +               ix8,iy8,iz8,ix9,iy9,iz9)
c
c     calculate the time at point 5 using stencil C and points 
c     0 through 4
c
      include 'fd.par'
      include 'fd.com'
c
      save=4./(vaverage(ix5,iy5,iz5,ix6,iy6,iz6)+
     +         vaverage(ix5,iy5,iz5,ix7,iy7,iz7)+
     +         vaverage(ix5,iy5,iz5,ix8,iy8,iz8)+
     +         vaverage(ix5,iy5,iz5,ix9,iy9,iz9))
c
      t0=time(ix0,iy0,iz0)
      t1=time(ix1,iy1,iz1)
      t2=time(ix2,iy2,iz2)
      t3=time(ix3,iy3,iz3)
      t4=time(ix4,iy4,iz4)
c
      t13=t1-t3
      t04=t0-t4
c
      sqrtterm=save*save-0.25*(t13*t13+t04*t04)
c
      if(sqrtterm.lt.0.) then          
        nnegsqrt=nnegsqrt+1
      else
        tc=t2+sqrt(sqrtterm)
        if(tc.lt.max(t0,t1,t3,t4)) then 
          ncaus=ncaus+1  
        else  
          if(tc.lt.tptmin) then
            tptmin=tc
            istencil=4
          end if
        end if
      end if
c
      nstencil=nstencil+1
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine timeb2d(ix1,iy1,iz1,ix2,iy2,iz2,
     +                 ix4,iy4,iz4,ix5,iy5,iz5,is)
c
c     calculate the time at point 5 using 2D version of stencil B and 
c     points 1, 2 and 4
c
      include 'fd.par'
      include 'fd.com'
c
      save=4./(vel(ix1,iy1,iz1)+vel(ix2,iy2,iz2)+
     +         vel(ix4,iy4,iz4)+vel(ix5,iy5,iz5))
c
      t1=time(ix1,iy1,iz1)
      t2=time(ix2,iy2,iz2)
      t4=time(ix4,iy4,iz4)
c
      t24=t2-t4
c
      sqrtterm=2.*save*save-t24*t24
c
      if(sqrtterm.lt.0.) then          
        nnegsqrt=nnegsqrt+1
      else
        tb2d=t1+sqrt(sqrtterm)
        if(tb2d.lt.max(t2,t4)) then
          ncaus=ncaus+1  
        else
          if(tb2d.lt.tptmin) then
            tptmin=tb2d
            istencil=is
          end if
        end if
      end if
c
      nstencil=nstencil+1
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine timec2d(ix1,iy1,iz1,ix2,iy2,iz2,ix4,iy4,iz4,
     +                 ix5,iy5,iz5,ix6,iy6,iz6,ix7,iy7,iz7)
c
c     calculate the time at point 5 using 2D version of stencil C and
c     points 1, 2 and 4
c
      include 'fd.par'
      include 'fd.com'
c
      save=8./(2.*(vel(ix5,iy5,iz5)+vel(ix2,iy2,iz2))+
     +             vel(ix6,iy6,iz6)+vel(ix1,iy1,iz1)+
     +             vel(ix7,iy7,iz7)+vel(ix4,iy4,iz4))
c
      t1=time(ix1,iy1,iz1)
      t2=time(ix2,iy2,iz2)
      t4=time(ix4,iy4,iz4)
c
      t14=t1-t4
c
      sqrtterm=save*save-0.25*t14*t14
c       
      if(sqrtterm.lt.0.) then          
        nnegsqrt=nnegsqrt+1
      else
        tc2d=t2+sqrt(sqrtterm)
        if(tc2d.lt.max(t1,t4)) then
          ncaus=ncaus+1  
        else
          if(tc2d.lt.tptmin) then
            tptmin=tc2d
            istencil=5
          end if
        end if
      end if
c
      nstencil=nstencil+1
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine time1d(ix1,iy1,iz1,ix2,iy2,iz2,is)
c
c     calculate the time at point 2 using 1-D stencil from  point 1
c
      include 'fd.par'
      include 'fd.com'
c
      save=2./(vel(ix1,iy1,iz1)+vel(ix2,iy2,iz2))
c
      t1d=time(ix1,iy1,iz1)+save
c
      if(t1d.lt.tptmin) then
        tptmin=t1d
        istencil=is
      end if
c
      nstencil=nstencil+1
c
      return
      end
