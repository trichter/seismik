c
c     version 1.0  Nov 1993
c
c     time routines for RAY
c
c     ----------------------------------------------------------------
c
      function timeinterp(xp,yp,zp,ixc,iyc,izc)
c
c     calculate the traveltime at the point (xp,yp,zp), in general not
c     a grid point
c
      include 'ray.par'
      include 'ray.com'
c        
      t1=time(ixc,iyc,izc)
      t2=time(ixc+1,iyc,izc)
      t5=time(ixc,iyc,izc+1)
      t6=time(ixc+1,iyc,izc+1)
c
      if(ny.eq.1) then
        t3=t1
        t4=t2
        t7=t5
        t8=t6
        yfrac=0.
      else
        t3=time(ixc,iyc+1,izc)
        t4=time(ixc+1,iyc+1,izc)
        t7=time(ixc,iyc+1,izc+1)
        t8=time(ixc+1,iyc+1,izc+1)
        yfrac=yp-(iyc-1)*size-ymin
      end if
c
      xfrac=xp-(ixc-1)*size-xmin
      zfrac=zp-(izc-1)*size-zmin
      yzfrac=yfrac*zfrac
c
      tleft= (size*((t3-t1)*yfrac+(t5-t1)*zfrac)+
     +       (t1-t3+t7-t5)*yzfrac)/size2+t1
      tright=(size*((t4-t2)*yfrac+(t6-t2)*zfrac)+
     +       (t2-t4+t8-t6)*yzfrac)/size2+t2
c
      timeinterp=(size-xfrac)*tleft+xfrac*tright
c
      return
      end
