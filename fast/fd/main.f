c
c     version 1.1  Apr 1995
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |            **************  F D  **************               |
c     |                                                              |
c     |     3D finite-difference traveltime modelling program        |
c     |                                                              |
c     |                 Written by Colin A. Zelt                     |
c     |                                                              |
c     |                Dept of Geology & Geophysics                  |   
c     |                      Rice University                         |
c     |                  Houston, TX  77005-1892                     | 
c     |                                                              |
c     ----------------------------------------------------------------
c
c	 patched for g77 compilation
c
c	 Drew Brenders (brenders@geoladm.geol.queensu.ca)
c 
c     I/O units:  
c                 
c        10 -- input:  program input parameters
c                 
c        11 -- output: summary iformation
c                 
c        17 -- output: calculated traveltime-distance pairs
c
c        19 -- output: all Calcomp plot calls
c
c        23 -- input:  velocity model in 3D binary format
c                 
c        27 -- input:  receiver locations at which to output calculated 
c                      times
c                 
c        28 -- output: receiver locations and calculated times
c                 
c        34 -- input:  parameters describing size of 3D grid
c                 
c        36 -- output: 3D traveltime field
c
c        37 -- input:  current iteration
c
c        38 -- input/output:  log file           
c                 
c
c     ----------------------------------------------------------------
c                 
c 
      program main
c
      include 'fd.par'
      include 'fd.com'                 
c
      integer isource(nsources),stencilt(8),reverse,rhold,pick
      integer*2 itime(nxmax)
      real xsource(nsources),ysource(nsources),zsource(nsources),
     +     hwmin,lambda,omin(nsources),omax(nsources),
     +     vabove(nsources),vbelow(nsources),ireal(nxmax)
c-kb
c      character*10 tfile,rfile,t2file
c      character*9 ofile
c      character*6 r2file
      character*15 tfile,rfile,t2file
      character*15 ofile
      character*15 r2file
c-kb
      character*1 alog
      character ifile*80
c                 
      namelist /pltpar/ iplot,imodxz,imodxy,imodyz,xtrms,
     +                  iout,inode,itimes,nodeht,igrid,
     +                  iroute,icol,iclear,iwrite,ttunc,ipar2,
     +                  iseg,xwndow,ywndow,colour,ibcol,ifcol
c
      namelist /axepar/ iaxlab,albht,xzxorig,xzyorig,
     +                  xyxorig,xyyorig,sep,
     +                  xtmin,xtmax,xmm,ntickx,ndecix,
     +                  ytmin,ytmax,ymm,nticky,ndeciy,
     +                  ztmin,ztmax,zmm,ntickz,ndeciz,
     +                  yzxorig,yzyorig
c                 
      namelist /propar/ i2d,tmax,istop,reverse,nreverse,iorder,
     +                  ttmin,tthwmin,hwmin,trmin,omin,omax,
     +                  interfce,ifile
c                 
      namelist /srcpar/ isource,xsource,ysource,zsource,nptsrc,inear,
     +                  vabove,vbelow,istype,pick
      namelist /spepar/ istype,iout,itimes,reverse,interfce
c
c     initialize parameters
c
c-kb
c      data xsource/nsources*-9999999./,
c     +     isource/nsources*-1/,tfile,rfile,r2file,ofile,t2file
c     +     /'fd  .times','fd  .picks','rec.  ','fd  .calc',
c     +      'fr  .times'/,omin/nsources*0./,omax/nsources*0./,
c     +     vabove/nsources*-1./,vbelow/nsources*-1./
c
      data xsource/nsources*-9999999./,
     +     isource/nsources*-1/,tfile,rfile,r2file,ofile,t2file
     +     /'fd    .times','fd    .picks','rec.    ','fd    .calc',
     +      'fr   .times'/,omin/nsources*0./,omax/nsources*0./,
     +     vabove/nsources*-1./,vbelow/nsources*-1./
c-kb
c
      write(6,335)
335   format(/'FD: finite difference traveltime calculation')
      open(50, file='nowrite', status='old', err=9999, iostat=ioflag)
9999  if(ioflag.gt.0) then
        io=0
      else
        io=51
      end if
c
      pick=1
      ipar2=0
      interfce=0
      xtrms=1000.
      istype=0
      nccrev=0
      trmin=.01
      hwmin=.01
      ttmin=.001
      tthwmin=.001
      iclear=0
      nreverse=10
      iorder=0
      reverse=7
      inear=0
      ttunc=.01
      istop=0
      iwrite=1
      tmax=100.
      do 1090 i=1,8
1090     stencilt(i)=0
      rmssum=0.
      chisum=0.
      nrms=0
      ntcalc=0
      ntplot=0
      ntstencil=0
      ntnegsqrt=0
      ntcaus=0
      nsrccalct=0
      itimes=0
      i2d=0
      inode=0
      isep=0
      iout=0
      iflagt=0
      nptsrc=5
      xzxorig=15.
      xzyorig=-1.
      xyxorig=-1.
      xyyorig=-1.
      yzxorig=-1.
      yzyorig=-1.
      iroute=1
      iaxlab=1
      imod=1
      imodxz=0
      imodxy=0
      imodyz=0
      nsource=0
c
      open(10, file='f.in', status='old')
c      
c     read in program control parameters from unit 10
c                 
      read(10,pltpar)
      read(10,axepar)
      read(10,propar)
      read(10,srcpar)
      close(10)
c
      if(ipar2.eq.1) then
        open(10, file='f2.in', status='old')
        read(10,spepar)
        close(10)
      end if
c
      if(iout.gt.1) then
        open(48, file='stop.in', status='old')
        read(48,*) iistop
        if(iistop.gt.0) stop
c
        open(47, file='time.out', status='old')
        read(47,*) iout
        rewind(47)
        write(47,*) -iout+1
        if(iout.eq.1) then
          open(46, file='lambda', status='old')
          read(46,*) lambda
          if(lambda.gt.0.) stop
        end if
      end if 
c
      if(i2d.ne.1.and.
     +  n2dmax.lt.max(nxmax*nymax,nxmax*nzmax,nymax*nzmax)) then
        write(0,1075)
1075    format(/'***  n2dmax set incorrectly  ***'/)
        stop
      end if
c
c     open I/O units
c
      if(iplot.eq.0) iplot=-1
      if(iplot.eq.2) iplot=0
      if(iplot.le.0) open(19, file='p.out')
c
      if(nptsrc.lt.3) then 
        write(0,795) 
795     format(/'***  nptsrc < 3  ***'/) 
        stop 
      end if 
c
      if(xzyorig.lt.0.) xzyorig=xzxorig
      if(yzxorig.lt.0.) yzxorig=xzxorig
      if(yzyorig.lt.0.) yzyorig=yzxorig
      if(xyxorig.lt.0.) xyxorig=xzxorig
      if(xyyorig.lt.0.) xyyorig=xyxorig
c
c     read in and setup the velocity model
c                 
      call model
c
      if(interfce.eq.1) then
        open(52, file=ifile, form='unformatted', status='old')
        do j=1,ny
           read(52) (ireal(i),i=1,nx)
           do i=1,nx
               inter(i,j)=nint(ireal(i))
           enddo
        enddo
      end if
c                 
c     calculate scale of each plot axis
c                 
      xscale=(xmax-xmin)/xmm
      yscale=(ymax-ymin)/ymm
      zscale=-(zmax-zmin)/zmm
      n2=nodeht/2.
      tmult=32766./tmax
      rhold=reverse
      n10=int(nodestotal/10.)
      if(reverse.eq.0) nreverse=0
      if(istop.eq.1) itimes=2
      if(reverse.ne.0.and.hwmin.gt.0.) then
        nhwmin=max(1,nint(nodestotal*hwmin/100.))
        write(io,875) nhwmin+1
875     format(
     +'minimum number of headwaves for reverse propagation:',i8)
      else
        nhwmin=0
      end if
c
c     set up the colours
c
      if(iroute.ne.1) then
        ibcol=0
      end if
      ncol=0
      do 1020 i=pcol,1,-1
         if(colour(i).ge.0) then
           ncol=i
           go to 1030
         end if
1020  continue
1030  if(ncol.eq.0) then
        colour(1)=2
        colour(2)=3
        colour(3)=4
        colour(4)=5
        colour(5)=6
        colour(6)=8
        colour(7)=17
        colour(8)=27
        colour(9)=22
        colour(10)=7
        ncol=10
      end if
      do 1040 i=1,8
         if(icol(i).lt.0) icol(i)=ifcol
1040  continue
c   
      if(i2d.eq.1) then
        imodxy=0
        imodyz=0
        if(iclear.eq.1) iclear=0
        if(abs(reverse).eq.3.or.abs(reverse).eq.4) then
          write(0,785)
785       format(/'***  invalid reverse side for 2D mode  ***'/)
          stop
        end if
      end if
      if(imodxz.eq.0.and.imodxy.eq.0.and.imodyz.eq.0) then
        noplot=1
        icflag=0
        iclear=0
      else
        noplot=0
      end if
c
      if(imodxz.eq.1) yzyorig=yzyorig+zmm+2.*sep
      if(imodxz.eq.1) xyyorig=xyyorig+zmm+sep
      if(imodyz.eq.1) xyyorig=xyyorig+zmm+2.*sep
c                 
c     determine nsource
c
      if(istype.ne.1) then
        do 1050 i=1,nsources
           if(xsource(i).le.-999999.) go to 1060
           nsource=nsource+1
           if(isource(i).gt.0.and.
     +       (xsource(i).lt.xmin.or.xsource(i).gt.xmax.or.
     +        ysource(i).lt.ymin.or.ysource(i).gt.ymax.or.
     +        zsource(i).lt.zmin.or.zsource(i).gt.zmax)) then
              write(0,4) i
4             format(/'***  source ',i3,' is outside model  ***'/)
                stop
         end if
1050    continue
      else
        do 1051 i=1,nsources
           if(isource(i).lt.0) go to 1060
           nsource=nsource+1
1051    continue
      end if
1060  if(nsource.eq.0) then
        write(0,1)
1       format(/'***  no sources specified  ***'/)
        stop
      end if
c
      if(inear.eq.1) then
        do i=1,nsource
          if(vabove(i).lt.0..or.vbelow(i).lt.0.) then
            write(0,515)
515         format(/'***  must specify vabove and vbelow  ***'/)
            stop
          end if
        enddo
      end if
c
      if(itimes.eq.1)     
     +open(27, file='rec.in', form='unformatted', status='old')
      if(itimes.eq.2) open(11, file='fd.out')
c
      ns=0
      do 1000 is=1,nsource
c
         if(isource(is).eq.0) go to 1000
c
         ns=ns+1
         if(iplots.eq.1) call aldone
c
         reverse=rhold
         nside=0
         nstop=0
         ncalc=0
         nsrccalc=0
         nplot=0
         nnegsqrt=0
         ncaus=0
         nstencil=0
         if(iclear.eq.0) then
           icflag=0
         else
           if(reverse.ne.0) then
             icflag=0
           else
             icflag=iclear
           end if
         end if
         do 1010 i=1,6
1010        tminside(i)=999999.
         do 1080 i=1,8
1080        stencil(i)=0
c         do 110 i=1,nxmax
c            do 110 k=1,nzmax
c110            xzplot(i,k)=0
c         do 140 i=1,nxmax
c            do 140 j=1,nymax
c140            xyplot(i,j)=0
c         do 120 j=1,nymax
c            do 120 k=1,nzmax
c120            yzplot(j,k)=0
c                 
c        plot velocity model
c
         if(imodxz.eq.1) call pltmodxz(igrid)
         if(imodyz.eq.1) call pltmodyz(imodxz,igrid)
         if(imodxy.eq.1) call pltmodxy(imodxz,imodyz,igrid)
c
         if(itimes.eq.2) then
c-kb
c           id1=is/10
c           id2=is-id1*10
c           rfile(3:4)=char(id1+48)//char(id2+48)
           id1=is/1000
           id2=(is-id1*1000)/100
           id3=(is-id1*1000-id2*100)/10
           id4=is-id1*1000-id2*100-id3*10   
           rfile(3:6)=char(id1+48)//char(id2+48)//char(id3+48)
     +//char(id4+48)
c-kb
c
           open(28, file=rfile, form='unformatted', status='old')
         end if
c
         if(istop.eq.1) then
           xsmin=xsource(is)
           xsmax=xsource(is)
           ysmin=ysource(is)
           ysmax=ysource(is)
           zsmin=zsource(is)
           zsmax=zsource(is)
180        read(28,end=900) xobs,yobs,zobs,tobs,uobs,inum
           if(inum.lt.0.or.inum.ne.pick) go to 180
           if(xobs.lt.xsmin) xsmin=xobs
           if(yobs.lt.ysmin) ysmin=yobs
           if(zobs.lt.zsmin) zsmin=zobs
           if(xobs.gt.xsmax) xsmax=xobs
           if(yobs.gt.ysmax) ysmax=yobs
           if(zobs.gt.zsmax) zsmax=zobs
           go to 180
900        rewind(28)
           call findnode(xsmin,ysmin,zsmin,ic1,jc1,kc1)
           call findnode(xsmax,ysmax,zsmax,ic2,jc2,kc2)
           sidestop(1)=max(1,ic1-1) 
           sidestop(2)=min(nx,ic2+1) 
           sidestop(3)=max(1,jc1-1) 
           sidestop(4)=min(ny,jc2+1) 
           sidestop(5)=max(1,kc1-1) 
           sidestop(6)=min(nz,kc2+1) 
         end if
c
c        initilize time field array to zero
c
         if(ny.ne.1) then
           do 10 i=0,nx+1
              do 10 j=0,ny+1
                 do 10 k=0,nz+1
10                  time(i,j,k)=1.e10
         else
           do 19 i=0,nx+1
                 do 19 k=0,nz+1
19                  time(i,1,k)=1.e10
         end if
c
         if(istype.ne.1) then
           if(noplot.eq.0) 
     +     call pltsource(xsource(is),ysource(is),zsource(is),
     +                    imodxz,imodxy,imodyz)
c
           call findnode(xsource(is),ysource(is),zsource(is),
     +                 ixsource,iysource,izsource)
c
           if(iwrite.eq.1)
     +     write(io,2) is,xsource(is),ysource(is),zsource(is),
     +              ixsource,iysource,izsource
2          format(/'source: ',i3,'  at (x,y,z): ',3f7.2/
     +           '          node location: ',3i7)
c
           if(noplot.eq.0)
     +     call pltsrcnode(ixsource,iysource,izsource,imodxz,imodxy,
     +                     imodyz)
c
           call timenearsrc(xsource(is),ysource(is),zsource(is),
     +          ixsource,iysource,izsource,nptsrc,imodxz,
     +          imodxy,imodyz,istop,inear,vabove(is),vbelow(is))
         else
           if(iwrite.eq.1)
     +     write(io,9) is
9          format(/'source: ',i3)
c
           call readsrc(is,imodxz,imodxy,imodyz,istop,tmult)
         end if
c
         if(istop.eq.1.and.iwrite.eq.1) write(io,905) 
     +     (sidestop(i),i=1,6)
905      format('          node propagation limit: ',6i5)
c
         if(iplots.eq.1) call empty
c
         if(iwrite.eq.1) write(io,8)
8        format(/
     +   '   % of nodes completed in forward propagation:'/
     +   '   -----------------------------------------------')
c
         if(i2d.ne.1) then
           call findiff(imodxz,imodxy,imodyz,istop,reverse,nreverse,
     +          iorder,ttmin,nhwmin,tthwmin,trmin,ncrev)
         else
           iys=int((ysource(is)-ymin)/size)+1
           ygrid=ymin+float(iys-1)*size
           if(abs(ygrid-ysource(is)).gt..001) then
             write(io,11)
11           format(/
     +       '***  source not on model grid in y direction  ***'/)
             go to 1000
           end if
c
           call findiff2d(imodxz,imodxy,imodyz,istop,iys,reverse,
     +     nreverse,iorder,ttmin,nhwmin,tthwmin,trmin,ncrev)
           do 1100 i=1,nx
              do 1100 j=1,ny
                 do 1100 k=1,nz
1100                time(i,j,k)=time(i,iys,k)
         end if
c
         if(nsrccalc.lt.ncalc) then
           if(iwrite.eq.1) write(io,6) 
     +     (stencil(i),i=1,8),
     +     (100.*float(stencil(i))/float(ncalc-nsrccalc),i=1,8)
6          format(/'number of nodes timed with each stencil:'/
     +     '     A: ',i7,'   B: ',i7,'  B2: ',i7,'   C: ',i7/
     +     '    C2: ',i7,'  C1: ',i7,'  D2: ',i7,'  D1: ',i7/
     +              'percentages:'/
     +     '     A: ',f7.3,'   B: ',f7.3,'  B2: ',f7.3,'   C: ',f7.3/
     +     '    C2: ',f7.3,'  C1: ',f7.3,'  D2: ',f7.3,'  D1: ',f7.3)
           do 1070 i=1,8
1070          stencilt(i)=stencilt(i)+stencil(i)
           nsrccalct=nsrccalct+nsrccalc
         else
           write(io,7)
7          format(/'all nodes timed within the source box')
         end if
c
         if(iwrite.eq.1) write(io,3) ncalc,nplot,nstencil,
     +                    nnegsqrt,ncaus,stencil(7)+stencil(8)
3        format('number of nodes timed:            ',i10/
     +          'number of nodes plotted:          ',i10/
     +          'number of stencils calculated:    ',i10/
     +          'number of negative square roots:  ',i10/
     +          'number of causality violations:   ',i10/
     +          'number of headwave stencils used: ',i10)
c
         if(reverse.ne.0) then
           write(io,1095) is,ncrev
1095       format('source ',i3,' completed using ',i3,' reversals')
           nccrev=nccrev+ncrev
         else
           write(io,1085) is
1085       format('source: ',i3,' completed')
         end if
c
         ntcalc=ntcalc+ncalc
         ntplot=ntplot+nplot
         ntstencil=ntstencil+nstencil
         ntnegsqrt=ntnegsqrt+nnegsqrt
         ntcaus=ntcaus+ncaus
c
         if(iout.eq.1.or.iout.eq.2.or.(iout.eq.3.and.istype.eq.0))
     +   then
c-kb
c           id1=is/10
c           id2=is-id1*10
c           tfile(3:4)=char(id1+48)//char(id2+48)
c           t2file(3:4)=char(id1+48)//char(id2+48)
           id1=is/1000
           id2=(is-id1*1000)/100
           id3=(is-id1*1000-id2*100)/10
           id4=is-id1*1000-id2*100-id3*10
           tfile(3:6)=char(id1+48)//char(id2+48)//char(id3+48)
     +//char(id4+48)
           t2file(3:6)=char(id1+48)//char(id2+48)//char(id3+48)
     +//char(id4+48)
c-kb
c
           if(istype.ne.1) then
             open(36, file=tfile, form='unformatted') 
           else
             open(36, file=t2file, form='unformatted') 
           end if
c
           do 131 k=1,nz
              do 131 j=1,ny
                 do 171 i=1,nx
                    if(time(i,j,k).gt.1.e9) time(i,j,k)=0.
171                 itime(i)=time(i,j,k)*tmult
131              write(36) (itime(i),i=1,nx)
           close(36)
         end if
c
         if(abs(itimes).eq.1) then
c-kb
c           id1=is/10
c           id2=is-id1*10
c           ofile(3:4)=char(id1+48)//char(id2+48)
           id1=is/1000
           id2=(is-id1*1000)/100
           id3=(is-id1*1000-id2*100)/10
           id4=is-id1*1000-id2*100-id3*10
           ofile(3:6)=char(id1+48)//char(id2+48)//char(id3+48)
     +//char(id4+48)
c-kb
c
           open(28, file=ofile, form='unformatted') 
c
           if(itimes.eq.-1) then
c-kb
c             r2file(5:6)=char(id1+48)//char(id2+48)
             r2file(5:8)=char(id1+48)//char(id2+48)//char(id3+48)
     +//char(id4+48)
c-kb
             open(27, file=r2file, form='unformatted') 
           else
             rewind(27)
           end if
150        read(27,end=1000) xrec,yrec,zrec,trec,urec,inum
           if(inum.eq.-2) then
             write(28) xrec,yrec,zrec,trec,urec,inum
             go to 1000
           end if
           if(inum.eq.-1) then
             write(28) xsource(is),ysource(is),zsource(is),0.,0.,inum
             go to 150
           end if
           if(pick.ne.inum) go to 150
           call find(xrec,yrec,zrec,irc,jrc,krc)
           if(ny.ne.1) then
           if(irc.lt.1.or.irc.gt.nx-1.or.jrc.lt.1.or.jrc.gt.ny-1.or.
     +        krc.lt.1.or.krc.gt.nz-1) go to 150
           trec=timeinterp(xrec,yrec,zrec,irc,jrc,krc)
           else  
           if(irc.lt.1.or.irc.gt.nx-1.or.
     +        krc.lt.1.or.krc.gt.nz-1) go to 150
           trec=timeinterp2d(xrec,zrec,irc,krc,1)
           end if
           if(itimes.eq.1) then
             ttuncr=ttunc
           else
             ttuncr=urec
           end if
           write(28) xrec,yrec,zrec,trec,ttuncr,inum
           go to 150
         end if
c
         if(itimes.eq.2) then
160        read(28,end=1000) xobs,yobs,zobs,tobs,uobs,inum
           if(inum.eq.-2) go to 1000
           if(inum.eq.-1) go to 160
           if(inum.ne.pick) go to 160
           if(omin(is).gt.0.or.omax(is).gt.0.) then
             dxy=((xsource(is)-xobs)**2+(ysource(is)-yobs)**2)**.5
             if(omin(is).gt.0..and.dxy.lt.omin(is)) go to 160
             if(omax(is).gt.0..and.dxy.gt.omax(is)) go to 160
           end if
c
           call find(xobs,yobs,zobs,irc,jrc,krc)
           if(ny.ne.1) then
           if(irc.lt.1.or.irc.gt.nx-1.or.jrc.lt.1.or.jrc.gt.ny-1.or.
     +        krc.lt.1.or.krc.gt.nz-1) go to 160
           tcalc=timeinterp(xobs,yobs,zobs,irc,jrc,krc)
           else  
           if(irc.lt.1.or.irc.gt.nx-1.or.
     +        krc.lt.1.or.krc.gt.nz-1) go to 160
           tcalc=timeinterp2d(xobs,zobs,irc,krc,1)
           end if
           diff2=(tobs-tcalc)*(tobs-tcalc)
           rmssum=rmssum+diff2
           chisum=chisum+diff2/uobs**2
           nrms=nrms+1
           go to 160
         end if
c
1000  continue
c
      open(37, file='current.iteration', status='old', err=99, 
     +     iostat=icistat)
c
99    if(icistat.gt.0) then
        open(37, file='current.iteration')
        ici=1
        write(37,*) ici
      else
        read(37,*) ici
      end if
c
      open(38, file='log.file')
155   read(38,55,end=98) alog
55    format(a1)
      go to 155
98    backspace(38) 
      write(38,65) ici
      write(io,65) ici
65    format(/'>>>  iteration: ',i5)
c
      write(38,33) ns,ntcalc,ntplot,ntstencil,ntnegsqrt,ntcaus,
     +             stencilt(7)+stencilt(8)
      if(rhold.ne.0) write(38,705) nccrev
705   format('total number of reverse propagations:   ',i10)
      if(ns.gt.0.and.ntstencil.gt.0) write(38,665)
     + (stencilt(i),i=1,8),
     + (100.*float(stencilt(i))/float(ntcalc-nsrccalct),i=1,8)
c
      if(ns.gt.1.and.ntstencil.gt.0) then
        write(io,33) ns,ntcalc,ntplot,ntstencil,ntnegsqrt,ntcaus,
     +              stencilt(7)+stencilt(8)
        if(rhold.ne.0) write(io,705) nccrev
33      format('------------------------------------------------'
     +,         '----------------------'/
     +       'total number of sources:                ',i10/
     +       'total number of nodes timed:            ',i10/
     +       'total number of nodes plotted:          ',i10/
     +       'total number of stencils calculated:    ',i10/
     +       'total number of negative square roots:  ',i10/
     +       'total number of causality violations:   ',i10/
     +       'total number of headwave stencils used: ',i10)
        if(nsrccalct.lt.ntcalc) then
          if(ns.gt.0.and.ntstencil.gt.0) write(io,665)
     +    (stencilt(i),i=1,8),
     +    (100.*float(stencilt(i))/float(ntcalc-nsrccalct),i=1,8)
665       format(
     + 'number of nodes timed with each stencil for all shots:'/
     +'     A: ',i8,'     B: ',i8,'    B2: ',i8,'     C: ',i8/
     +'    C2: ',i8,'    C1: ',i8,'    D2: ',i8,'    D1: ',i8/
     + 'percentages:'/
     +'     A: ',f7.3,'     B: ',f7.3,'    B2: ',f7.3,'     C: ',f7.3/
     +'    C2: ',f7.3,'    C1: ',f7.3,'    D2: ',f7.3,'    D1: ',f7.3)
        end if
      end if
      write(io,695)
695   format(/)
c
      if(nrms.gt.1) then
        trms=xtrms*(rmssum/float(nrms))**.5
        chi=chisum/float(nrms-1)
        write(io,125) nrms,trms,chi
        write(38,125) nrms,trms,chi
        write(11,125) nrms,trms,chi
125     format('RMS traveltime residual for ',i8,'  data: ',f12.5/
     +         '                    Normalized chi-squared: ',f12.4)
      end if
c
      if(iplots.eq.1) then
        call empty
        call plotnd(1)
      end if
c
      stop        
      end         
