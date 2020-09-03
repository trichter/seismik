c
c     version 1.1  Apr 1995
c
c     2D Finite difference routine for FD
c
c     ----------------------------------------------------------------
c
      subroutine findiff2d(imodxz,imodxy,imodyz,istop,j,reverse,
     +           nreverse,iorder,ttmin,nhwmin,tthwmin,trmin,ncrev)
c
c     time all remaining points outside the source box using 2D finite
c     difference solution of eikonal equation
c
      include 'fd.par'
      include 'fd.com'
c
      integer ix(n1dmax),iz(n1dmax),ipos(n1dmax),reverse,
     +        revopp(6),hface(6),hpface(6)
      real times(n1dmax)
      character ans*1,face(6)*6
      data revopp/2,1,4,3,6,5/
     +     face/'  left',' right','  back',' front','   top','bottom'/
      ncrev=0
      ireverse=0
      trmax=1.e10
      ioc=0
      do 1010 i=1,6
1010     hwside(i)=0
      hface(1)=nx
      hface(2)=1
      hface(5)=nz
      hface(6)=1
c
      if(side(3).ne.1) then
        side(3)=1
        nside=nside+1
        nstop=nstop+1
      end if
      if(side(4).ne.ny) then
        side(4)=ny
        nside=nside+1
        nstop=nstop+1
      end if
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
        z1=float(side(5)-1)*size+zmin
        z2=float(side(6)-1)*size+zmin
        xp1=(x1-xmin)/xscale+xzxorig
        xp2=(x2-xmin)/xscale+xzxorig
        zp1=(z1-zmax)/zscale+xzyorig
        zp2=(z2-zmax)/zscale+xzyorig
        call box(xp1,zp1,xp2,zp2)
        call empty
      end if
c
      call expand(iorder,istop,ioc,tmincurrent,nt,nrepmem,nreplace,
     +            nfrev)
c
      go to (1111,2222,9998,9998,5555,6666) iside
c      
9998  write(6,95)
95    format(/'***  could not decide which side to expand  ***'/)
      stop
c
1111  continue
c
c       expand left side of box
c
        ixp=side(1)
        do 10 k=side(5),side(6)
           nt=nt+1
           iz(nt)=k
10         times(nt)=time(ixp,j,k)
c
        call sort(times,ipos,nt)
c
        do 110 i=1,nt
c
1000       iyp=j
           izp=iz(ipos(i))
           tptmin=time(ixp-1,iyp,izp)
           torig=tptmin
           ipltnode=5
           istencil=0
c
           call stencils2d(ixp,iyp,izp,ttrans)
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
        do 20 k=side(5),side(6)
           nt=nt+1
           iz(nt)=k
20         times(nt)=time(ixp,j,k)
c
        call sort(times,ipos,nt)
c
        do 220 i=1,nt
c
2000       iyp=j 
           izp=iz(ipos(i))
           tptmin=time(ixp+1,iyp,izp)
           torig=tptmin
           ipltnode=5
           istencil=0
c
           call stencils2d(ixp,iyp,izp,ttrans)
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
5555  continue
c
c       expand top side of box
c
        izp=side(5)
        do 50 i=side(1),side(2)
           nt=nt+1
           ix(nt)=i
50         times(nt)=time(i,j,izp)
c
        call sort(times,ipos,nt)
c
        do 550 i=1,nt
c
5000       ixp=ix(ipos(i))
           iyp=j
           tptmin=time(ixp,iyp,izp-1)
           torig=tptmin
           ipltnode=5
           istencil=0
c
           call stencils2d(ixp,iyp,izp,ttrans)
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
           nt=nt+1
           ix(nt)=i
60         times(nt)=time(i,j,izp)
c
        call sort(times,ipos,nt)
c
        do 660 i=1,nt
c
6000       ixp=ix(ipos(i))
           iyp=j
           tptmin=time(ixp,iyp,izp+1)
           torig=tptmin
           ipltnode=5
           istencil=0
c
           call stencils2d(ixp,iyp,izp,ttrans)
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
c
      call pass(nreplace,treplace,trmax,ncrev,reverse,nreverse,
     +          nhwmin,trmin,istop,hface,hpface,nfrev,irev)   
c
      if(irev.eq.1) go to 10000
c
      return
      end
