c 
c     version 1.1  Apr 1995
c
c     3D Finite difference routine for FD
c
c     ---------------------------------------------------------------- 
c
      subroutine findiff(imodxz,imodxy,imodyz,istop,reverse,nreverse,
     +           iorder,ttmin,nhwmin,tthwmin,trmin,ncrev)
c
c     time all remaining points outside the source box using finite
c     difference solution of eikonal equation
c
      include 'fd.par'
      include 'fd.com'
c
      integer ix(n2dmax),iy(n2dmax),iz(n2dmax),ipos(n2dmax),reverse,
     +        hface(6),hpface(6)
      real times(n2dmax)
      character ans*1
      ncrev=0
      ireverse=0
      trmax=1.e10
      ioc=0
      do 1010 i=1,6
1010     hwside(i)=0
      hface(1)=nx
      hface(2)=1
      hface(3)=ny
      hface(4)=1
      hface(5)=nz
      hface(6)=1
c
10000 continue
c
      if(nside.eq.6.or.(istop.eq.1.and.nstop.eq.6)) then
        if(ireverse.eq.1) then
           if(iwrite.eq.1) write(io,785) side(iside)
785        format('     last face:     ',i5)
        end if
        go to 9999
      end if
c
      if(icflag.gt.1) then
        write(6,*) 'iside=',iside,'   side=',side(iside)
        read(5,1) ans
1       format(a1)
        call erase
        call pcolor(icol(3))
        x1=float(side(1)-1)*size+xmin
        x2=float(side(2)-1)*size+xmin
        y1=float(side(3)-1)*size+ymin
        y2=float(side(4)-1)*size+ymin
        z1=float(side(5)-1)*size+zmin
        z2=float(side(6)-1)*size+zmin
        xp1=(x1-xmin)/xscale+xzxorig
        xp2=(x2-xmin)/xscale+xzxorig
        zp1=(z1-zmax)/zscale+xzyorig
        zp2=(z2-zmax)/zscale+xzyorig
        if(imodxz.eq.1) call box(xp1,zp1,xp2,zp2)
        yp1=(y1-ymin)/yscale+yzxorig
        yp2=(y2-ymin)/yscale+yzxorig
        zp1=(z1-zmax)/zscale+yzyorig
        zp2=(z2-zmax)/zscale+yzyorig
        if(imodyz.eq.1) call box(yp1,zp1,yp2,zp2)
        xp1=(x1-xmin)/xscale+xyxorig
        xp2=(x2-xmin)/xscale+xyxorig
        yp1=-(y1-ymax)/yscale+xyyorig
        yp2=-(y2-ymax)/yscale+xyyorig
        if(imodxy.eq.1) call box(xp1,yp1,xp2,yp2)
        call empty
      end if
c
      call expand(iorder,istop,ioc,tmincurrent,nt,nrepmem,nreplace,
     +            nfrev)
c
      go to (1111,2222,3333,4444,5555,6666) iside
c      
      write(6,95)
95    format(/'***  could not decide which side to expand  ***'/)
      stop   
c
1111  continue
c
c       expand left side of box
c
        ixp=side(1)
        do 10 j=side(3),side(4)
           do 10 k=side(5),side(6)
              nt=nt+1
              iy(nt)=j
              iz(nt)=k
10            times(nt)=time(ixp,j,k)
c
        call sort(times,ipos,nt)
c
        do 110 i=1,nt
c
           iyp=iy(ipos(i))
           izp=iz(ipos(i))
           if(interfce.eq.1.and.izp.gt.inter(ixp-1,iyp)) go to 110
           tptmin=time(ixp-1,iyp,izp)
           torig=tptmin               
           ipltnode=5
           istencil=0
c
           call stencils(ixp,iyp,izp,ttrans)
c
           call fnode(ixp-1,iyp,izp,torig,ttrans,ttmin,trmax,
     +  tthwmin,nreplace,treplace,hface,imodxz,imodxy,imodyz)
110     continue
c
        call fside(nfrev,nreplace,nrepmem,-1,ixp,hpface(iside),irep)
c
        if(irep.eq.0) go to 9999
        go to 10000
c      
2222  continue
c
c       expand right side of box
c
        ixp=side(2)
        do 20 j=side(3),side(4)
           do 20 k=side(5),side(6)
              nt=nt+1
              iy(nt)=j
              iz(nt)=k
20            times(nt)=time(ixp,j,k)
c
        call sort(times,ipos,nt)
c
        do 220 i=1,nt
c
           iyp=iy(ipos(i))
           izp=iz(ipos(i))
           if(interfce.eq.1.and.izp.gt.inter(ixp+1,iyp)) go to 220
           tptmin=time(ixp+1,iyp,izp)
           torig=tptmin               
           ipltnode=5
           istencil=0
c
           call stencils(ixp,iyp,izp,ttrans)
c
           call fnode(ixp+1,iyp,izp,torig,ttrans,ttmin,trmax,
     +  tthwmin,nreplace,treplace,hface,imodxz,imodxy,imodyz)
220     continue
c
        call fside(nfrev,nreplace,nrepmem,1,-ixp,-hpface(iside),irep)
c
        if(irep.eq.0) go to 9999
        go to 10000
c      
3333  continue
c
c       expand back side of box
c
        iyp=side(3)
        do 30 i=side(1),side(2)
           do 30 k=side(5),side(6)
              nt=nt+1
              ix(nt)=i
              iz(nt)=k
30            times(nt)=time(i,iyp,k)
c
        call sort(times,ipos,nt)
c
        do 330 i=1,nt
c
           ixp=ix(ipos(i))
           izp=iz(ipos(i))
           if(interfce.eq.1.and.izp.gt.inter(ixp,iyp-1)) go to 330
           tptmin=time(ixp,iyp-1,izp)
           torig=tptmin               
           ipltnode=5
           istencil=0
c
           call stencils(ixp,iyp,izp,ttrans)
c
           call fnode(ixp,iyp-1,izp,torig,ttrans,ttmin,trmax,
     +  tthwmin,nreplace,treplace,hface,imodxz,imodxy,imodyz)
330     continue
c
        call fside(nfrev,nreplace,nrepmem,-1,iyp,hpface(iside),irep)
c
        if(irep.eq.0) go to 9999
        go to 10000
c      
4444  continue
c
c       expand front side of box
c
        iyp=side(4)
        do 40 i=side(1),side(2)
           do 40 k=side(5),side(6)
              nt=nt+1
              ix(nt)=i
              iz(nt)=k
40            times(nt)=time(i,iyp,k)
c
        call sort(times,ipos,nt)
c
        do 440 i=1,nt
c
           ixp=ix(ipos(i))
           izp=iz(ipos(i))
           if(interfce.eq.1.and.izp.gt.inter(ixp,iyp+1)) go to 440
           tptmin=time(ixp,iyp+1,izp)
           torig=tptmin               
           ipltnode=5
           istencil=0
c
           call stencils(ixp,iyp,izp,ttrans)
c
           call fnode(ixp,iyp+1,izp,torig,ttrans,ttmin,trmax,
     +  tthwmin,nreplace,treplace,hface,imodxz,imodxy,imodyz)
440     continue
c
        call fside(nfrev,nreplace,nrepmem,1,-iyp,-hpface(iside),irep)
c
        if(irep.eq.0) go to 9999
        go to 10000
c      
5555  continue
c
c       expand top side of box
c
        izp=side(5)
        do 50 i=side(1),side(2)
           do 50 j=side(3),side(4)
              nt=nt+1
              ix(nt)=i
              iy(nt)=j
50            times(nt)=time(i,j,izp)
c
        call sort(times,ipos,nt)
c
        do 550 i=1,nt
c
           ixp=ix(ipos(i))
           iyp=iy(ipos(i))
           if(interfce.eq.1.and.izp-1.gt.inter(ixp,iyp)) go to 550
           tptmin=time(ixp,iyp,izp-1)
           torig=tptmin               
           ipltnode=5
           istencil=0
c
           call stencils(ixp,iyp,izp,ttrans)
c
           call fnode(ixp,iyp,izp-1,torig,ttrans,ttmin,trmax,
     +  tthwmin,nreplace,treplace,hface,imodxz,imodxy,imodyz)
550     continue
c
        call fside(nfrev,nreplace,nrepmem,-1,izp,hpface(iside),irep)
c
        if(irep.eq.0) go to 9999
        go to 10000
c      
6666  continue
c
c       expand bottom side of box
c
        izp=side(6)
        do 60 i=side(1),side(2)
           do 60 j=side(3),side(4)
              nt=nt+1
              ix(nt)=i
              iy(nt)=j
60            times(nt)=time(i,j,izp)
c
        call sort(times,ipos,nt)
c
        do 660 i=1,nt
c
           ixp=ix(ipos(i))
           iyp=iy(ipos(i))
           if(interfce.eq.1.and.izp+1.gt.inter(ixp,iyp)) go to 660
           tptmin=time(ixp,iyp,izp+1)
           torig=tptmin               
           ipltnode=5
           istencil=0
c
           call stencils(ixp,iyp,izp,ttrans)
c
           call fnode(ixp,iyp,izp+1,torig,ttrans,ttmin,trmax,
     +  tthwmin,nreplace,treplace,hface,imodxz,imodxy,imodyz)
660     continue
c
        call fside(nfrev,nreplace,nrepmem,1,-izp,-hpface(iside),irep)
c
        if(irep.eq.0) go to 9999
        go to 10000
c
9999  continue
c
      call pass(nreplace,treplace,trmax,ncrev,reverse,nreverse,
     +          nhwmin,trmin,istop,hface,hpface,nfrev,irev)
c
      if(irev.eq.1) go to 10000
c
      return
      end
