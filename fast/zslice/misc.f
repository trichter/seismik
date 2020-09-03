c
c     version 1.0  Feb 1995
c
c     misc routines for ZSLICE
c
c     ----------------------------------------------------------------
c
      function timeinterp(xp,yp,zp,ixc,iyc,izc)
c
c     calculate the traveltime at the point (xp,yp,zp), in general not
c     a grid point
c
      include 'zslice.par'
      include 'zslice.com'
c
      t1=time(ixc,iyc,izc)
      t2=time(ixc+1,iyc,izc)
      t3=time(ixc,iyc+1,izc)
      t4=time(ixc+1,iyc+1,izc)
      t5=time(ixc,iyc,izc+1)
      t6=time(ixc+1,iyc,izc+1)
      t7=time(ixc,iyc+1,izc+1)
      t8=time(ixc+1,iyc+1,izc+1)
c
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
      if(t1.eq.undef.or.t2.eq.undef.or.t3.eq.undef.or.t4.eq.undef.or.
     +   t5.eq.undef.or.t6.eq.undef.or.t7.eq.undef.or.t8.eq.undef)
     +   timeinterp=undef
c
      return
      end
