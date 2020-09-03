c
c     version 1.1  Apr 1995
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |            **************  R A Y  *************              |
c     |                                                              |
c     |       Calculate and plot ray paths given a 3D time grid      |
c     |       and calculate slowness and interfce perturbations     |
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
c        16 -- output: receiver locations and calculated traveltimes
c                 
c        17 -- input/output: slowness perturbation sums
c                 
c        18 -- input/output: weight (inverse pick uncertainty) sums
c                 
c        19 -- output: all Calcomp plot calls
c
c        20 -- input: receiver locations and observed traveltimes
c
c        29 -- input/output: ray counts in each model cell
c
c        35 -- input:  parameters describing size of 3D grid
c
c        36 -- input:  3D time grid from the program FD
c
c        37 -- input:  2D interfce above which velocities are fixed
c
c        38 -- input/output:  log file
c
c        45 -- output:  chi**2 value
c                 
c
c     ----------------------------------------------------------------
c                 
c 
      program main
c
      include 'ray.par'
      include 'ray.com'                 
c 
      integer source(6),itrace(nsources),npti(2),s4,oophist(100),
     +        moophist(100),k2(nxmax,nymax),pick,outray
      real xpl(2),ypl(2),lambda,p1(3),p2(3),p3(3),p12(3),p13(3),
     +     l1oop(100),l2oop(100),msoop,nix,niy,niz
c-kb
c      character ierr*1,reply*1,tfile*10,rfile*10,rtfile*11,tlab*12,
c     +          intfile*80,r2file*9,t2file*10
      character ierr*1,reply*1,tfile*15,rfile*15,rtfile*15,tlab*15,
     +          intfile*80,r2file*15,t2file*15
c-kb
c
      namelist /pltpar/ iplot,iroute,iray,ircol,iscol,itrms,irec,irefr,
     +                  iseg,xwndow,ywndow,ibcol,ifcol,colour,i2d,
     +                  symht,igrid,sep,istep,iwarn,iout,itomo,souht,
     +                  ireccol,ndecim,ndecir,ires,xres,yres,zres,
     +                  iscreen,ixy,ixz,iyz,i3d,theta,npskip,imid,hwo,
     +                  iwater,ifirst,xtrms,imethod,iunc,uave,ioop,
     +                  noopmin,ipar2,outray
c
      namelist /axepar/ iaxlab,albht,xorig,yorig,xomm,
     +                  xmm,ymm,zmm,tmm,xomin,xomax,ttmin,ttmax,
     +                  xotmin,xotmax,tttmin,tttmax,ndecixo,ndecit,
     +                  ntickxo,ntickt
c
      namelist /raypar/ nptsrc,smin,itrace,tmax,omin,omax,interfce,
     +                  nptmax,istype,intfile,pick
c
      namelist /ttpar/ itime,vred,itccol,itocol,itrcol,ttunc
      namelist /spepar/ istype,itrms
c
c     initialize parameters
c
c-kb
c      data tfile,rfile,rtfile,r2file,t2file/'fd  .times','fd  .picks',
c     +     'ray  .times','fd  .refl','fr  .times'/,
c     +      itrace/nsources*0/,oophist/100*0/,moophist/100*0/
      data tfile,rfile,rtfile,r2file,t2file/'fd    .times',
     +     'fd    .picks',
     +     'ray    .times','fd    .refl','fr    .times'/,
     +      itrace/nsources*0/,oophist/100*0/,moophist/100*0/
c-kb
c
      write(6,335)
335   format('RAY: ray tracing')
      open(50, file='nowrite', status='old', err=9999, iostat=ioflag)
9999  if(ioflag.gt.0) then
        io=0
      else
        io=51
      end if
c
      pick=1
      ipar2=0
      istype=0
      ioop=0
      noopmin=10
      iyz=0
c M.P. 18.04.2011 change maximum number of node points
c                 which a ray can cross
      outray=0
      nptmax=nxmax+nymax+2*nzmax
      uave=.05
      nlc=0
      nk=0
      nl=0
      xtrms=1000.
      iwater=0
      irefr=0
      ifirst=0
      hwo=60.
      imid=0
      npskip=1
      theta=20.
      i3d=0
      ttunc=.01
      ixy=1
      ixz=1
      interfce=0
      iscreen=1
      i2d=0
      vred=0.
      itccol=-1
      itocol=-1
      itrcol=-1
      ires=0
      omin=-1.
      omax=-1.
      ndecim=1
      ndecir=1
      ireccol=1
      irec=0
      ncount=0
      nrayshot=0
      ntrayshot=0
      itrms=0
      rmssum=0.
      chisum=0.
      rmssm2=0.
      chism2=0.
      num2=0
      rmssm3=0.
      chism3=0.
      num3=0
      nsource=0
      ntfail=0
      ntptst=0
      nttrace=0
      ntpmax=0
      itomo=0
      iout=0
      iwarn=0
      smin=.001
      iscol=1
      nptsrc=10
      istep=0
      ircol=1
      iray=1
      igrid=0
      imethod=0
      iunc=0
      symht=.1
      souht=1.
      isep=0
      tmax=100.
      xorig=15.
      yorig=-1.
      iroute=1
      iaxlab=1
      xomin=-999999.
      xomax=-999999.
      ttmin=-999999.
      ttmax=-999999.
      noop=0
      soop=0.
      msoop=0.
c
      open(10, file='r.in', status='old')
c      
c     read in program control parameters from unit 10
c                 
      read(10,pltpar)
      read(10,axepar)
      read(10,raypar)
      read(10,ttpar)
      close(10)
c
      if(ipar2.eq.1) then
        open(10, file='r2.in', status='old')
        read(10,spepar)
        close(10)
      end if
c
      if(itomo.ge.5) then
        open(48, file='stop.in', status='old')
        read(48,*) iistop
        if(iistop.gt.0) stop
c
        open(46, file='lambda', status='old')
        read(46,*) lambda
        if(lambda.gt.0.) stop
      end if
c
      if(itomo.eq.2) then
        open(27, file='tomo', status='old')
        read(27,*) itomo
        if(itomo.ne.1.and.itomo.ne.2) then
          write(0,675)
675       format(/'***  invalid value for itomo  ***'/)
          stop
        end if
        rewind(27)   
        write(27,*) -itomo+3
        close(27)
      end if
c
      if(ires.eq.1.and.itomo.eq.1) itomo=3
      if(abs(itime).gt.0) iray=0
      if(i3d.eq.1) then
        theta=theta/57.29577951
        cost=cos(theta)
        sint=sin(theta)
        ixy=0
        ixz=1
      end if
c
c     open I/O units
c
      open(15, file='r.out')
      if(iplot.eq.0) iplot=-1
      if(iplot.eq.2) iplot=0
      if(iplot.le.0) open(19, file='p.out')
      open(35, file='for.header', status='old')
      if(itrms.gt.0) open(45, file='chi')
c
      if(nptsrc.lt.3) then
        write(0,795)
795     format(/'***  nptsrc < 3  ***'/)
        stop
      end if
c
      read(35,115) xmin,xmax,ymin,ymax,zmin,zmax,size,nx,ny,nz
115   format(7f10.3,3i10)
c
      if(float(nx**2+ny**2+nz**2)**.5.gt.nray) write(0,755)
755   format(/'***  nray should be greater than model diagonal  ***'/)
c
      nodestotal=nx*ny*nz
      write(io,25) xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,nodestotal
25    format(/'model dimensions:'/
     +        '-----------------'/
     +        'xmin, xmax: ',2f10.3/'ymin, ymax: ',2f10.3/,
     +        'zmin, zmax: ',2f10.3/
     +        'number of model nodes in each direction (x,y,z):',3i5/
     +        'total number of model nodes:                    ',i15)
      write(io,35) size
35    format('node spacing: ',f7.3,' km')
c
      if(nx.gt.nxmax.or.ny.gt.nymax.or.nz.gt.nzmax) then
        write(0,3)
3       format(/'***  model is too large  ***'/)
        stop
      end if
c
      if(nptmax.gt.0) 
     +  open(49, file='rays.xyz', status='unknown')
c     +  open(49, file='rays.xyz', form='unformatted')
c
      if(itomo.gt.0) then
        open(39, file='inv.header', status='old')
        read(39,*) inx,iny,inz
        xii=(xmax-xmin)/float(inx)
        yii=(ymax-ymin)/float(iny)
        zii=(zmax-zmin)/float(inz)
        nm=inx*iny*inz
      end if
c
      if( itomo.ge.5) then
        if(istype.ne.2) 
     +    open(40, file='delta.time', form='unformatted')
        if(istype.eq.0) then
          open(41, file='data.kernel', form='unformatted')
          open(42, file='row.kernel', form='unformatted')
          open(43, file='column.kernel', form='unformatted')
          open(44, file='nzero.kernel')
        else
          if(istype.eq.1) then
            open(41, file='data1.kernel', form='unformatted')
            open(42, file='row1.kernel', form='unformatted')
            open(43, file='column1.kernel', form='unformatted')
            open(44, file='nzero1.kernel')
          else
            if(itomo.eq.5) then
              open(41, file='data2.kernel', form='unformatted')
              open(42, file='row2.kernel', form='unformatted')
              open(43, file='column2.kernel', form='unformatted')
              open(44, file='nzero2.kernel')
            else
              open(52, file='int.deriv', form='unformatted')
            end if
          end if
        end if
        open(29,file='num.cell',form='unformatted',status='old')
        do 2131 k=1,inz
           do 2131 j=1,iny
2131          read(29) (numray(i,j,k),i=1,inx)
      end if
c
      if(interfce.eq.1) then
        open(37, file='bathymetry', form='unformatted',
     +         status='old')
        do 1300 j=1,ny
           read(37) (boundary(i,j),i=1,nx)
           do 1310 i=1,nx
1310          boundary(i,j)=boundary(i,j)/1000.
1300    continue
      end if
c
      if(istype.gt.0) then
        open(54, file=intfile, form='unformatted', status='old')
           do j=1,ny
              read(54) (inter(i,j),i=1,nx)
              do i=1,nx
                  k2(i,j)=nint(inter(i,j))
              enddo
           enddo
           do j=1,ny
              read(54) (inter(i,j),i=1,nx)
              do i=1,nx
                  inter(i,j)=inter(i,j)/1000.
              enddo
           enddo
      end if
c
      if(iray.eq.0) istep=0
      if(itrms.eq.2) iray=0
c
c     calculate scale of each plot axis
c
      xscale=(xmax-xmin)/xmm
      yscale=-(ymax-ymin)/ymm
      zscale=-(zmax-zmin)/zmm
      if(yorig.lt.0.) yorig=xorig
      xzxorig=xorig
      xzyorig=yorig
      xyxorig=xzxorig
      xyyorig=xzyorig+zmm+sep
      sizex4=size*4.
      sizex2=size*2.
      size2=size*size
      sized2=size*.5
      tmult=32766./tmax
      txsize=tmult*size
      txsizex2=tmult*sizex2
      one=1
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
c
      do 10 i=1,nx
10       xmodel(i)=xmin+float(i-1)*size
      do 20 j=1,ny
20       ymodel(j)=ymin+float(j-1)*size
      do 30 k=1,nz
30       zmodel(k)=zmin+float(k-1)*size
c
      if(itomo.eq.1) then
        open(17,file='sl.sums', form='unformatted', status='old')
        open(18,file='weight.cell',form='unformatted',status='old')
        open(29,file='num.cell',form='unformatted',status='old')
        if(imethod.eq.1) then
        open(30, file='slow.mod', form='unformatted', status='old')
c
          do 310 k=1,inz
             do 310 j=1,iny
310              read(30) (imodel(i,j,k),i=1,inx)
        end if
c
        do 2130 k=1,inz
           do 2130 j=1,iny
              read(18) (pert2(i,j,k),i=1,inx)
              read(29) (numray(i,j,k),i=1,inx)
2130          read(17) (pert1(i,j,k),i=1,inx)
      end if
c
      if(itomo.eq.2) then
        open(17,file='sl.pert', form='unformatted', status='old')
        open(28,file='delta.times', form='unformatted')
        do 2140 k=1,inz
           do 2140 j=1,iny
2140          read(17) (pert1(i,j,k),i=1,inx)
      end if
c
      if(itomo.eq.3) then
        open(17, file='res.ker', form='unformatted', status='old')
        do 2170 k=1,inz
           do 2170 j=1,iny
2170          read(17) (pert1(i,j,k),i=1,inx)
        call findnode(xres,yres,zres,ixres,iyres,izres)
        ixres1=ixres-1
        iyres1=iyres-1
        izres1=izres-1
      end if
c
      do 4000 iss=1,nsources
c
      nfail=0
      nrayshot=0
      nptst=0
      ntrace=0
      npmax=0
c
      if(itrace(iss).lt.1.or.itrace(iss).gt.nsources) go to 4000
      is=itrace(iss) 
c
      nsource=nsource+1
c-kb
c      id1=is/10
c      id2=is-id1*10
c      rfile(3:4)=char(id1+48)//char(id2+48)
c      if(irec.ge.2) go to 7000
c      tfile(3:4)=char(id1+48)//char(id2+48)
c      t2file(3:4)=char(id1+48)//char(id2+48)
c      rtfile(4:5)=char(id1+48)//char(id2+48)
c      r2file(3:4)=char(id1+48)//char(id2+48)
      id1=is/1000
      id2=(is-id1*1000)/100
      id3=(is-id1*1000-id2*100)/10
      id4=is-id1*1000-id2*100-id3*10
      rfile(3:6)=char(id1+48)//char(id2+48)//char(id3+48)
     +//char(id4+48)
      if(irec.ge.2) go to 7000
      tfile(3:6)=char(id1+48)//char(id2+48)//char(id3+48)
     +//char(id4+48)
      t2file(3:6)=char(id1+48)//char(id2+48)//char(id3+48)
     +//char(id4+48)
      rtfile(4:7)=char(id1+48)//char(id2+48)//char(id3+48)
     +//char(id4+48)
      r2file(3:6)=char(id1+48)//char(id2+48)//char(id3+48)
     +//char(id4+48)
c-kb
      if(iout.eq.1) open(16, file=rtfile, form='unformatted')
      if(istype.eq.1) open(53, file=r2file, form='unformatted')
c
      write(15,95)
95    format(/'source    xrec     yrec     zrec   i   j   k'
     +,'      x        y        z   npts prob')
c
      if(istype.eq.1) then
        open(36, file=t2file, form='unformatted', status='old')
      else
        open(36, file=tfile, form='unformatted', status='old')
      end if
      do 130 k=1,nz
         do 130 j=1,ny
130         read(36) (time(i,j,k),i=1,nx)
c
      if(istype.eq.1) then
        do i=1,nx
           do j=1,ny
              if(k2(i,j).lt.nz) then
                 idt=time(i,j,k2(i,j)-1)-time(i,j,k2(i,j))
                 do k=k2(i,j)+1,nz
                    time(i,j,k)=time(i,j,k-1)-idt
                 enddo
              end if
           enddo
         enddo
      end if
c
      close(36)
c
7000  if(iray.gt.0) then
        if(iplots.eq.0) then
          call plots(xwndow,ywndow,iroute)
          call segmnt(1)
          iplots=1
        else
          if(iray.eq.1) go to 5000
          call empty
          call aldone
          call erase
        end if
c
        call pcolor(ifcol)
        if(ixz.eq.1.and.iyz.ne.1) 
     +    call box(xzxorig,xzyorig,xzxorig+xmm,xzyorig+zmm)
        if(ixz.eq.1.and.iyz.eq.1) 
     +    call box(xzxorig,xzyorig,xzxorig+ymm,xzyorig+zmm)
        if(i3d.eq.1) then
          xo=xzxorig+ymm*cost
          zo=xzyorig+ymm*sint
          call box(xo,zo,xo+xmm,zo+zmm)
          call plot(xzxorig,xzyorig,3)
          call plot(xo,zo,2)
          call plot(xzxorig+xmm,xzyorig,3)
          call plot(xo+xmm,zo,2)
          call plot(xzxorig,xzyorig+zmm,3)
          call plot(xo,zo+zmm,2)
          call plot(xzxorig+xmm,xzyorig+zmm,3)
          call plot(xo+xmm,zo+zmm,2)
        end if
        if(ixy.eq.1) call box(xyxorig,xyyorig,xyxorig+xmm,xyyorig+ymm)
c
        if(iaxlab.eq.1) then
          call axtick(xmin,xmax,xtmin,xtmax,ntickx,ndecix)
          call axtick(zmin,zmax,ztmin,ztmax,ntickz,ndeciz)
          call axtick(ymin,ymax,ytmin,ytmax,nticky,ndeciy)
c
          if(ixz.eq.1) then
            if(iyz.ne.1) then
            call axis(xzxorig,xzyorig,xmin,xmax,xmm,xscale,0.,1,
     +           xtmin,xtmax,ntickx,ndecix,'X (km)',6,albht)
            else
            call axis(xzxorig,xzyorig,ymin,ymax,ymm,-yscale,0.,1,
     +           ytmin,ytmax,nticky,ndeciy,'Y (km)',6,albht)
            end if
          else
            call axis(xyxorig,xyyorig,xmin,xmax,xmm,xscale,0.,1,
     +           xtmin,xtmax,ntickx,ndecix,'X (km)',6,albht)
          end if
          if(ixz.eq.1)
     +       call axis(xzxorig,xzyorig,zmin,zmax,zmm,zscale,90.,1,
     +       ztmin,ztmax,ntickz,ndeciz,'Z (km)',6,albht)
          if(ixy.eq.1)
     +       call axis(xyxorig,xyyorig,ymin,ymax,ymm,yscale,90.,1,
     +       ytmin,ytmax,nticky,ndeciy,'Y (km)',6,albht)
c         call plotxy
        end if
c
        if(igrid.eq.1) then
          if(ixz.eq.1) then
            do 60 i=1,nx,ndecim
               xpl(1)=size*float(i-1)/xscale+xzxorig
               xpl(2)=xpl(1)
               ypl(1)=xzyorig
               ypl(2)=xzyorig+zmm
60             call line(xpl,ypl,2)
            do 70 k=1,nz,ndecim
             ypl(1)=(zmin+size*float(k-1)-zmax)/zscale+xzyorig
               ypl(2)=ypl(1)
               xpl(1)=xzxorig
               xpl(2)=xzxorig+xmm
70             call line(xpl,ypl,2)
          end if
          if(ixy.eq.1) then
            do 80 i=2,nx-1,ndecim
               xpl(1)=size*float(i-1)/xscale+xyxorig
               xpl(2)=xpl(1)
               ypl(1)=xyyorig
               ypl(2)=xyyorig+ymm
80             call line(xpl,ypl,2)
            do 90 j=2,ny-1,ndecim
               ypl(1)=(ymin+size*float(j-1)-ymax)/yscale+xyyorig
               ypl(2)=ypl(1)
               xpl(1)=xyxorig
               xpl(2)=xyxorig+xmm
90             call line(xpl,ypl,2)
          end if
        end if
c
        call empty
      end if
c
      if(abs(itime).gt.0) then
        if(iplots.eq.0) then
          call plots(xwndow,ywndow,iroute)
          call segmnt(1)
          iplots=1
          if(xomin.lt.-999998.) xomin=xmin
          if(xomax.lt.-999998.) xomax=xmax
          if(ttmin.lt.-999998.) ttmin=-tmax
          if(ttmax.lt.-999998.) ttmax=tmax
          tscale=(ttmax-ttmin)/tmm
          xoscale=(xomax-xomin)/xomm
          xtxorig=xorig
          xtyorig=yorig
          if(vred.le.0) then
            rvred=0.
          else
            rvred=1./vred
          end if
          if(itccol.lt.0) itccol=ifcol+1
          if(itocol.lt.0) itocol=ifcol
          if(itrcol.lt.0) itrcol=ifcol
        else
          if(itime.gt.0) go to 5000
          call empty
          call aldone
          call segmnt(1)
          call erase
        end if
c
        call pcolor(ifcol)
        call box(xtxorig,xtyorig,xtxorig+xomm,xtyorig+tmm)
        if(abs(itime).eq.2) then
          tt0=-ttmin/tscale+xtyorig
          call plot(xtxorig,tt0,3)
          call plot(xtxorig+xomm,tt0,2)
        end if
c
        if(iaxlab.eq.1) then
          call axtick(xomin,xomax,xotmin,xotmax,ntickxo,ndecixo)
          call axis(xtxorig,xtyorig,xomin,xomax,xomm,xoscale,0.,1,
     +       xotmin,xotmax,ntickxo,ndecixo,'offset (km)',11,albht)
          call axtick(ttmin,ttmax,tttmin,tttmax,ntickt,ndecit)
          if(abs(itime).eq.1) then
            if(vred.eq.0.) then
              nchart=8
              tlab='time (s)'
            else
              i1=int(vred)
              id1=int((vred-float(i1))*10+.05)
              id2=int((vred-float(i1)-float(id1)/10.)*100.+.5)
              if(id1.eq.0.and.id2.eq.0) then
                nchart=9
                tlab='T-D/  (s)'
                tlab(5:5)=char(i1+48)
              else
                if(id2.eq.0) then
                  nchart=11
                  tlab='T-D/    (s)'
                  tlab(5:7)=char(i1+48)//'.'//char(id1+48)
                else
                  nchart=12
                  tlab='T-D/     (s)'
                  tlab(5:8)=char(i1+48)//'.'//char(id1+48)//char(id2+48)
                end if
              end if
            end if
            call axis(xtxorig,xtyorig,ttmin,ttmax,tmm,tscale,90.,1,
     +      tttmin,tttmax,ntickt,ndecit,tlab,nchart,albht)
          else
            call axis(xtxorig,xtyorig,ttmin,ttmax,tmm,tscale,90.,1,
     +      tttmin,tttmax,ntickt,ndecit,'traveltime residual (s)',23,
     +      albht)
          end if
        end if
      end if
c
5000  if(istype.ne.2) then
        open(20, file=rfile, form='unformatted', status='old')
      else
        open(20, file=r2file, form='unformatted', status='old')
      end if
c
1000  read(20,end=999) xr,yr,zr,tr,ur,icol
c
      if(icol.eq.-2) go to 999
      iplotr=0
      if(icol.eq.-1.and.iout.eq.1) write(16) xr,yr,zr,0.,0.,-1
      if(icol.eq.-1.and.istype.eq.1) then
        write(53) xr,yr,zr,0.,0.,-1
        if(iray.gt.0.and.ircol.ge.0) call pcolor(ircol)
        go to 1000
      end if
      if(icol.eq.-1.and.istype.ne.1) then
c
        xs=xr
        ys=yr
        zs=zr
c
c       find closest node to the source
c
        call findnode(xs,ys,zs,ixs,iys,izs)
c
        if(ixs.lt.1.or.ixs.gt.nx.or.iys.lt.1.or.iys.gt.ny.or.
     +     izs.lt.1.or.izs.gt.nz) then
          write(io,895) is
895       format(/'***  source ',i3,' is outside model  ***'/)
          isrc=0
          go to 1000
        end if
        isrc=1
c
        il=(nptsrc-1)/2
c
        if(ixs-il.lt.1) then
          source(1)=1
        else
          source(1)=ixs-il
        end if
        if(ixs+il.gt.nx) then
          source(2)=nx
        else
          source(2)=ixs+il
        end if
        if(iys-il.lt.1) then
          source(3)=1
        else
          source(3)=iys-il
        end if
        if(iys+il.gt.ny) then
          source(4)=ny
        else
          source(4)=iys+il
        end if
        if(izs-il.lt.1) then
          source(5)=1
        else
          source(5)=izs-il
        end if
        if(izs+il.gt.nz) then
          source(6)=nz
        else
          source(6)=izs+il
        end if
c
c       plot source box
c
        if(iray.gt.0.and.iscol.ge.0.and.iscol.ne.ibcol) 
     +    call pltsrcbox(source,iscol,ixy,ixz,i3d)
c
        if(iray.gt.0.and.ircol.ge.0) call pcolor(ircol)
c
c       find which model cell the source is in
c
        call find(xs,ys,zs,ixs,iys,izs)
c
        go to 1000

      end if
      if(icol.ne.pick) go to 1000
c
      if(omin.gt.0.or.omax.gt.0.) then
        offset=((xr-xs)**2+(yr-ys)**2)**.5
        if(omax.gt.0..and.offset.gt.omax) go to 1000
        if(omin.gt.0..and.offset.lt.omin) go to 1000
      end if
c
      if(irec.ge.2) then
        ntrace=ntrace+1
        if(mod(ntrace-1,ndecir).ne.0) go to 1000
        if(irec.eq.2.and.nsource.ne.1.and.ifirst.eq.1) go to 1000
        xplot(1)=(xr-xmin)/xscale+xzxorig
        yplot(1)=(yr-ymax)/yscale+xyyorig
        zplot(1)=(zr-zmax)/zscale+xzyorig
        if(irec.eq.2.and.ireccol.ne.ibcol) then
          if(ixz.eq.1) call dot(xplot(1),zplot(1),symht,ireccol)
          if(ixy.eq.1) call dot(xplot(1),yplot(1),symht,ireccol)
        else
          xplot(2)=(xs-xmin)/xscale+xzxorig
          yplot(2)=(ys-ymax)/yscale+xyyorig
          zplot(2)=(zs-zmax)/zscale+xzyorig
          if(imid.le.0) then
c           if(ixz.le.0) call line(xplot,zplot,2)
            if(ixz.eq.1) call line(xplot,zplot,2)
            if(ixy.eq.1) call line(xplot,yplot,2)
          else
            if(imid.eq.1) then
              xplot(1)=(xplot(1)+xplot(2))/2. 
              yplot(1)=(yplot(1)+yplot(2))/2. 
              if(ixy.eq.1) call dot(xplot(1),yplot(1),symht,ircol)
            else
              dist=((xr-xs)**2+(yr-ys)**2)**.5
              xplot(1)=(xr+(xs-xr)*hwo/dist-xmin)/xscale+xzxorig
              yplot(1)=(yr+(ys-yr)*hwo/dist-ymax)/yscale+xyyorig
              xplot(2)=(xs+(xr-xs)*hwo/dist-xmin)/xscale+xzxorig
              yplot(2)=(ys+(yr-ys)*hwo/dist-ymax)/yscale+xyyorig
              if(ixy.eq.1) call line(xplot,yplot,2)
            end if
          end if
        end if 
        go to 1000
      end if
c
      if(istype.ne.1.and.isrc.eq.0) go to 1000
c
      npts=1
      xray(1)=xr
      yray(1)=yr
      zray(1)=zr
      iok=1
      ierr=' '
      xp=xray(npts)
      yp=yray(npts)
      zp=zray(npts)
c
c     find which model cell the receiver is in
c
      call find(xp,yp,zp,ic,jc,kc)
c
      cell(npts,1)=ic
      cell(npts,2)=jc
      cell(npts,3)=kc
c
c     check to see if receiver is outside the model
c
      if(ic.lt.1.or.ic.gt.nx-1.or.jc.lt.1.or.
     +  (jc.gt.ny-1.and.ny.gt.1).or.kc.lt.1.or.kc.gt.nz-1) then
        iok=0
        ierr='r'
        go to 2000
      end if
c
c     check to see if receiver is at the source
c
      if(istype.ne.1) then
        if(abs(xray(npts)-xs).lt.smin.and.abs(yray(npts)-ys).
     +    lt.smin.and.abs(zray(npts)-zs).lt.smin) then
          iok=0
          ierr='0'
          go to 2000
        end if
      end if
c
      if(istype.eq.2) then
        x1=xmin+(ic-1)*size
        x2=x1+size
        y1=ymin+(jc-1)*size
        y2=y1+size
        z1=inter(ic,jc)
        z2=inter(ic+1,jc)
        z3=inter(ic,jc+1)
        z4=inter(ic+1,jc+1)
        c1=y2*(z2-z1)-y1*(z4-z3)
        c3=x1*z2-x2*z1+x2*z3-x1*z4
        c4=z1-z2+z4-z3
        c7=size2
        a1=c1/c7
        a2=c3/c7
        a3=c4/c7

        call grad(xr,yr,zr,ar,br,cr)
        nix=a1+a3*yr
        niy=a2+a3*xr
        niz=-1.
        aray=sqrt(ar*ar+br*br+cr*cr)
        anorm=sqrt(nix*nix+niy*niy+niz*niz)
        costhe=(nix*ar+niy*br+niz*cr)/(aray*anorm)
        cosalp=1./anorm
        write(52) xr,yr,zr,tr,ur,ic,jc,costhe,cosalp
c       write(59,375) xr,yr,zr,tr,ur,ic,jc,costhe,cosalp
375     format(5f9.3,2i6,2f9.3)
      end if
c
      if(itrms.eq.2) then
c        print*,xr,yr,zr,cell(1,1),cell(1,2),cell(1,3)
        tcalc=timeinterp(xr,yr,zr,cell(1,1),cell(1,2),cell(1,3))
     +        /txsize
        diff2=(tr-tcalc)*(tr-tcalc)
        rmssum=rmssum+diff2
        chisum=chisum+diff2/ur**2
        if(abs(itime).gt.0.and.mod(ntrace-1,ndecir).eq.0) then
c         offset=((xr-xs)**2+(yr-ys)**2+(zr-zs)**2)**.5
          offset=((xr-xs)**2+(yr-ys)**2)**.5
          xop=(offset-xomin)/xoscale+xtxorig
          if(abs(itime).eq.1) then
            ttpo=(tr-offset*rvred-ttmin)/tscale+xtyorig
            ttpc=(tcalc-offset*rvred-ttmin)/tscale+xtyorig
            call dot(xop,ttpo,symht,itocol)
            call dot(xop,ttpc,symht,itccol)
          else
            ttpr=(tcalc-tr-ttmin)/tscale+xtyorig
            call dot(xop,ttpr,symht,itrcol)
          end if
        end if
        ntrace=ntrace+1
        go to 1000
      end if
c
c     check to see if receiver is within source box
c
      if(istype.ne.1) then
        if(ic.ge.source(1).and.ic.lt.source(2).and.
     +     jc.ge.source(3).and.jc.lt.source(4).and.
     +     kc.ge.source(5).and.kc.lt.source(6)) then
          iok=-1
          go to 2000
        end if
      end if
c
3000  continue
c
c     calculate normal vector to local traveltime field 
c     for current ray point
c
c      print*,"GRAD: ",xp,yp,zp,xg,yg,zg
      call grad(xp,yp,zp,xg,yg,zg)
c
c     calculate intersection point of ray with cell boundary
c
c      print*,"INT1: ",ic,jc,kc,xs,ys,zs,ixs,iys,izs
c      print*,"INT1: ",xg,yg,zg,xp,yp,zp,xi,yi,zi
      call intersect(ic,jc,kc,xs,ys,zs,ixs,iys,izs,xg,yg,zg,
     +               xp,yp,zp,xi,yi,zi,1)
c
c      print*,"INT2: ",ic,jc,kc,xs,ys,zs,ixs,iys,izs
c      print*,"INT2: ",xg,yg,zg,xp,yp,zp,xi,yi,zi
      npts=npts+1
      xray(npts)=xi
      yray(npts)=yi
      zray(npts)=zi
      xp=xi
      yp=yi
      zp=zi
c
      cell(npts,1)=ic
      cell(npts,2)=jc
      cell(npts,3)=kc
c
c     check to see if ray is within source box
c
      if(istype.ne.1) then
        if(ny.ne.1) then
          s4=source(4)
        else
          s4=2
        end if
        if(ic.ge.source(1).and.ic.lt.source(2).and.
     +     jc.ge.source(3).and.jc.lt.s4.and.
     +     kc.ge.source(5).and.kc.lt.source(6)) go to 2000
      else
        icr=cell(npts-1,1)
        jcr=cell(npts-1,2)
        kcr=cell(npts-1,3)
        zb=bndinterp(xp,yp,icr,jcr,2)
        if(zb.le.zp) then    
          x1=xmin+(icr-1)*size
          x2=x1+size
          y1=ymin+(jcr-1)*size
          y2=y1+size
          z1=inter(icr,jcr)  
          z2=inter(icr+1,jcr)  
          z3=inter(icr,jcr+1)  
          z4=inter(icr+1,jcr+1)  
          c1=y2*(z2-z1)-y1*(z4-z3)
          c3=x1*z2-x2*z1+x2*z3-x1*z4
          c4=z1-z2+z4-z3
          c5=y2*(x2*z1-x1*z2)-y1*(x2*z3-x1*z4)
          c7=size2
          a1=c1/c7
          a2=c3/c7
          a3=c4/c7
          a4=c5/c7
          ar=xray(npts)-xray(npts-1)
          br=yray(npts)-yray(npts-1)
          cr=zray(npts)-zray(npts-1)
          x0=xray(npts-1)
          y0=yray(npts-1)
          z0=zray(npts-1)
          aq=a3*ar*br
          bq=-cr+a1*ar+a2*br+a3*(ar*y0+br*x0)
          cq=-z0+a1*x0+a2*y0+a3*x0*y0+a4
          if(aq.eq.0.) then
            t=-cq/bq
            xi=x0+t*ar 
            yi=y0+t*br 
            zi=z0+t*cr 
          else
            det=max(0.,bq**2-4.*aq*cq)
            sdet=sqrt(det)
            t1=(-bq-sdet)/(2.*aq)
            t2=(-bq+sdet)/(2.*aq)
            if(abs(t1).le.abs(t2)) then
              xi=x0+t1*ar 
              yi=y0+t1*br 
              zi=z0+t1*cr 
            else
              xi=x0+t2*ar 
              yi=y0+t2*br 
              zi=z0+t2*cr 
            end if
          end if
          xray(npts)=xi
          yray(npts)=yi
          zray(npts)=zi
          go to 2000
        end if
      end if
c
c     check to see if ray has hit the edge of the model
c
      if(xp.le.xmin.or.xp.ge.xmax.or.(yp.le.ymin.and.ny.gt.1).or.
     +(yp.ge.ymax.and.ny.gt.1).or.zp.le.zmin.or.zp.ge.zmax) then
c
c       ray has reached the edge of the model
c
        iok=0
        ierr='>'
        go to 2000
      end if
c
      if(npts.eq.nray) then
c
c       ray consists of too many points
c
2400    iok=0
        ierr='!'
        go to 2000
      end if
c
      go to 3000
c
2000  continue
c
c     ray has either reached source box (iok=1) or has failed to do 
c     so (iok=0) or the ray started within the source box (iok=-1)
c
      if(abs(iok).eq.1) then
c
        if(istype.ne.1) then
c
c         calculate straight-line ray path inside source box to the 
c         source
c
          call straight(ic,jc,kc,xs,ys,zs,ixs,iys,izs,iflag)
          if(iflag.eq.1) goto 2400
        end if
c
        if(iout.eq.1.or.itomo.gt.0.or.itrms.eq.1.or.abs(itime).gt.0) 
     +  then
c          print*,xr,yr,zr,cell(1,1),cell(1,2),cell(1,3)
          tcalc=timeinterp(xr,yr,zr,cell(1,1),cell(1,2),cell(1,3))
     +          /txsize
          if(istype.eq.1.and.zi.gt.zmin.and.zi.lt.zmax) 
     +        write(53) xi,yi,zi,tr-tcalc,ur,icol
          if(iout.eq.1) write(16) xr,yr,zr,tcalc,ttunc,icol
          if(itrms.eq.1) then
            diff2=(tr-tcalc)*(tr-tcalc)
            rmssum=rmssum+diff2
            chisum=chisum+diff2/ur**2
          end if
          if(abs(itime).gt.0.and.mod(ntrace-1,ndecir).eq.0) then
c           offset=((xr-xs)**2+(yr-ys)**2+(zr-zs)**2)**.5
            offset=((xr-xs)**2+(yr-ys)**2)**.5
            xop=(offset-xomin)/xoscale+xtxorig
            if(abs(itime).eq.1) then
              ttpo=(tr-offset*rvred-ttmin)/tscale+xtyorig
              ttpc=(tcalc-offset*rvred-ttmin)/tscale+xtyorig
              call dot(xop,ttpo,symht,itocol)
              call dot(xop,ttpc,symht,itccol)
            else
              ttpr=(tcalc-tr-ttmin)/tscale+xtyorig
              call dot(xop,ttpr,symht,itrcol)
            end if
            write(34,395) offset,tr-tcalc,ur,(tr-tcalc)/ur
395         format(4f10.3)
          end if
c
          if(itomo.eq.1.or.itomo.eq.2.or.itomo.ge.5) then
            if(itomo.eq.1)
     +        call backproj(tr,ur,tcalc,interfce,iflag)
            if(itomo.eq.2)
     +        call d_time(tr,ur,tcalc,interfce,iflag)
            if(itomo.eq.5) 
     +       call kernel(tr,ur,uave,tcalc,interfce,iflag)
            if(itrms.eq.1.and.iflag.eq.0) then
              num2=num2+1
              rmssm2=rmssm2+diff2
              chism2=chism2+diff2/ur**2
            end if
          else
            if(itomo.eq.3) 
     +        call resolution(ixres,iyres,izres,ixres1,iyres1,izres1)
          end if
        end if
      else
        if(ierr.ne.'0') then
          nfail=nfail+1
        else
          nrayshot=nrayshot+1
        end if
        if(iwarn.eq.1) 
     +    write(io,85) is,xr,yr,zr,ic,jc,kc,xray(npts),
     +                yray(npts),zray(npts),npts,ierr 
        write(15,85) is,xr,yr,zr,ic,jc,kc,xray(npts),
     +               yray(npts),zray(npts),npts,ierr 
85      format(i5,3f9.3,3i4,3f9.3,i5,1x,a1)
      end if
c
2010  ntrace=ntrace+1
      nptst=nptst+npts
      if(npts.gt.npmax.and.npts.ne.nray) npmax=npts
c
      if(irefr.eq.1.and.iplotr.ne.1) go to 1200
      if(iok.gt.0.and.iray.gt.0.and.mod(ntrace-1,ndecir).eq.0) then
c     if(iok.ge.0.and.iray.gt.0.and.mod(ntrace-1,ndecir).eq.0) then
        nplt=0
        if(irefr.eq.1) then
          np1=npti(1)
          np2=npti(2)+1
        else
          np1=1
          np2=npts
        end if
        if(i3d.ne.1) then
          do 160 i=np1,np2-1,npskip
             nplt=nplt+1
             xplot(nplt)=(xray(i)-xmin)/xscale+xzxorig
             yplot(nplt)=(yray(i)-ymax)/yscale+xyyorig
             zplot(nplt)=(zray(i)-zmax)/zscale+xzyorig
160       continue
          nplt=nplt+1
          xplot(nplt)=(xray(np2)-xmin)/xscale+xzxorig
          yplot(nplt)=(yray(np2)-ymax)/yscale+xyyorig
          zplot(nplt)=(zray(np2)-zmax)/zscale+xzyorig
        else
          do 162 i=np1,np2-1,npskip
             nplt=nplt+1
             yplot(nplt)=(yray(i)-ymax)/yscale
             xplot(nplt)=(xray(i)-xmin)/xscale+xzxorig+yplot(nplt)*cost
             zplot(nplt)=(zray(i)-zmax)/zscale+xzyorig+yplot(nplt)*sint
162       continue
          nplt=nplt+1
          yplot(nplt)=(yray(np2)-ymax)/yscale
          xplot(nplt)=(xray(np2)-xmin)/xscale+
     +                 xzxorig+yplot(nplt)*cost
          zplot(nplt)=(zray(np2)-zmax)/zscale+
     +                 xzyorig+yplot(nplt)*sint
        end if
c
c       plot the ray path
c
c     write(47,*) xray(1),-zray(1)
c     write(47,*) xray(npts),-zray(npts)
c     write(47,*) '>'
c
        if(ibcol.ne.ircol) then
          if(nplt.gt.1) then
            if(ircol.lt.0) call pcolor(icol)
            if(ixz.eq.1.and.iyz.ne.1) call line(xplot,zplot,nplt)
            if(ixy.eq.1) call line(xplot,yplot,nplt)
            if(ixz.eq.1.and.iyz.eq.1) then
              nplt=0
              do 161 i=np1,np2-1,npskip
                 nplt=nplt+1
                 yplot(nplt)=-(yray(i)-ymin)/yscale+xzxorig
161           continue
              nplt=nplt+1
              xplot(nplt)=(xray(np2)-xmin)/xscale+xzxorig
              yplot(nplt)=-(yray(np2)-ymin)/yscale+xzxorig
              call line(yplot,zplot,nplt)
            end if
          end if
        end if
c
        if(irec.gt.0.and.ireccol.ne.ibcol) then
          if(ixz.eq.1) call dot(xplot(1),zplot(1),symht,ireccol)
          if(ixy.eq.1) call dot(xplot(1),yplot(1),symht,ireccol)
          if(ircol.ge.0) call pcolor(ircol)
        end if
      end if
c
      if(ioop.eq.1.and.npts.ge.noopmin) then 
        p1(1)=xray(1)
        p1(2)=yray(1)
        p1(3)=zray(1)
        p2(1)=xray(npts)
        p2(2)=yray(npts)
        p2(3)=zray(npts)
        p3(1)=(p1(1)+p2(1))/2.
        p3(2)=(p1(2)+p2(2))/2.
        p3(3)=14.
        p12(1)=p2(1)-p1(1)
        p12(2)=p2(2)-p1(2)
        p12(3)=p2(3)-p1(3)
        p13(1)=p3(1)-p1(1)
        p13(2)=p3(2)-p1(2)
        p13(3)=p3(3)-p1(3)
        a=p12(2)*p13(3)-p12(3)*p13(2)
        b=p12(3)*p13(1)-p12(1)*p13(3)
        c=p12(1)*p13(2)-p12(2)*p13(1)
        d=-a*p1(1)-b*p1(2)-c*p1(3)
        rl=sqrt((p1(1)-p2(1))**2+(p1(2)-p2(2))**2+(p1(3)-p2(3))**2)
        denom=sqrt(a**2+b**2+c**2)
        sum=0.
        dmax=0.
        do i=2,npts-1
           doff=abs(a*xray(i)+b*yray(i)+c*zray(i)+d)
           if(doff.gt.dmax) dmax=doff
           sum=sum+doff
        enddo
        sum=sum/denom
        dmax=dmax/denom
        ave=sum/(npts-2)
        rave=ave/rl
        rdmax=dmax/rl
        soop=soop+rave
        noop=noop+1
        msoop=msoop+rdmax
        nrave=nint(100.*rave)
        nrdmax=nint(100.*rdmax)
        oophist(nrave+1)=oophist(nrave+1)+1
        moophist(nrdmax+1)=moophist(nrdmax+1)+1
      end if
c
1200  if(istep.eq.1.and.mod(ntrace-1,ndecir).eq.0) then
        write(0,65) ntrace,xr,yr,zr
65      format('ray: ',i8,'   receiver: ',3f9.3,' ',$)
        call empty
        read(5,15) reply
15      format(a1)
        if(reply(1:1).eq.'s') go to 999
        if(reply(1:1).eq.'0') then
          istep=0
        else
          if(npts.gt.1) then
            call pcolor(ibcol)
            if(ixz.eq.1.and.iyz.ne.1) call line(xplot,zplot,npts)
            if(ixz.eq.1.and.iyz.eq.1) then
              do 1095 i=1,npts
1095             yplot(i)=-(yray(i)-ymin)/yscale+xzxorig
              call line(yplot,zplot,npts)
            end if
            if(ixy.eq.1.and.iyz.eq.1) then
              do 1096 i=1,npts
1096             yplot(i)=(yray(i)-ymax)/yscale+xyyorig
            end if
            if(ixy.eq.1) call line(xplot,yplot,npts)
            call pcolor(ircol)
          end if
        end if
      end if
      if(nptmax.gt.0 .and. outray.eq.1) then
        write(49,*)iss,xs,ys,zs,ntrace,xr,yr,zr
        do i=1,npts
         write(49,*) xray(i),yray(i),zray(i)
        enddo
      end if

      go to 1000
c
999   continue 
c
      if(istype.ne.1.and.iray.gt.0.and.iscol.le.0.and.
     +  iscol.ne.ibcol) then
        xplot(1)=(xs-xmin)/xscale+xzxorig
        yplot(1)=(ys-ymax)/yscale+xyyorig
        zplot(1)=(zs-zmax)/zscale+xzyorig
        if(ixy.eq.1) call dot(xplot(1),yplot(1),souht,abs(iscol))
        if(ixz.eq.1.and.iyz.ne.1) then
          call dot(xplot(1),zplot(1),souht,abs(iscol))
        else
          if(iyz.eq.1) then 
            yplot(1)=-(ys-ymin)/yscale+xzxorig
            call dot(yplot(1),zplot(1),souht,abs(iscol))
          end if
        end if
      end if
c
      if(ntrace.gt.0) then
        npave=nint(float(nptst)/float(ntrace))
      else
        npave=0
      end if
c
      if(iscreen.eq.1) then
        write(io,75) is,ntrace,nptst,nfail,nrayshot,npave,npmax
      else
        write(io,1085) is,ntrace
1085    format('source: ',i3,' completed tracing ',i8,' rays')
      end if
c
      write(15,75) is,ntrace,nptst,nfail,nrayshot,npave,npmax
75    format(/
     +  'source number:                                 ',i8/
     +  'number of rays/points traced:                  ',i8,i10/
     +  'number of rays failed to reach source:         ',i8/
     +  'number of rays starting at source location:    ',i8/
     +  'average/maximum number of points defining ray: ',i8,i10)
c
      ntrayshot=ntrayshot+nrayshot
      ntfail=ntfail+nfail
      ntptst=ntptst+nptst
      nttrace=nttrace+ntrace
      if(npmax.gt.ntpmax) ntpmax=npmax
c
4000  continue
c
      if(itomo.ne.2) then
        open(38, file='log.file')
155     read(38,55,end=98) alog
55      format(a1)
        go to 155
      end if
c
98    if(nttrace.gt.0) then
        npave=nint(float(ntptst)/float(nttrace))
      else
        npave=0
      end if
c
      if(itomo.ne.2)then
       backspace(38)
       write(38,135) nsource,nttrace,ntptst,ntfail,
     +               ntrayshot,npave,ntpmax
      endif
135     format(/
     +  'total number of sources:                          ',i8/
     +  'total number of rays/points traced:               ',i8,i10/
     +  'total number of rays failed to reach source:      ',i8/
     +  'total number of rays starting at source location: ',i8/
     +  'average/maximum number of points defining ray:    ',i8,i10/)
c
      if(nsource.gt.1) then
        write(io,105) nsource,nttrace,ntptst,ntfail,ntrayshot,npave,
     +               ntpmax
        write(15,105) nsource,nttrace,ntptst,ntfail,ntrayshot,npave,
     +                ntpmax
105     format(/
     +  '---------------------------------------------------------',
     +  '---------------'/
     +  'total number of sources:                          ',i8/
     +  'total number of rays/points traced:               ',i8,i10/
     +  'total number of rays failed to reach source:      ',i8/
     +  'total number of rays starting at source location: ',i8/
     +  'average/maximum number of points defining ray:    ',i8,i10/)
      end if
c
      if(itrms.gt.0.and.nttrace.gt.ntfail+ntrayshot+1) then
        trms=xtrms*(rmssum/float(nttrace-ntfail-ntrayshot))**.5
        chi=chisum/float(nttrace-ntfail-ntrayshot-1)
        if(nsource.eq.1) then
          write(io,127) 
          write(15,127)
          if(itomo.ne.2) write(38,127)
127       format(' ')
        end if
        write(io,125) nttrace-ntfail-ntrayshot,trms,chi
        write(15,125) nttrace-ntfail-ntrayshot,trms,chi
        if(itomo.ne.2) 
     +    write(38,125) nttrace-ntfail-ntrayshot,trms,chi
125     format('RMS traveltime residual for ',i8,' rays: ',f12.5/
     +       '                   Normalized chi-squared: ',f12.4/)
        if(num2.le.1.or.iwater.ne.1) write(45,*) chi
        if(num2.gt.1.and.iwater.eq.1) then
          trms2=xtrms*(rmssm2/float(num2))**.5
          chi2=chism2/float(num2-1)
          write(io,129) num2,trms2,chi2
          write(15,129) num2,trms2,chi2
          if(itomo.ne.2) write(38,129) num2,trms2,chi2
129       format('excluding water wave arrivals:'/
     +           'RMS traveltime residual for ',i8,' rays: ',f12.5/
     +         '                   Normalized chi-squared: ',f12.4/)
          write(45,*) chi2
        end if
        if(num3.gt.1) then
          trms3=xtrms*(rmssm3/float(num3))**.5
          chi3=chism3/float(num3-1)
          write(io,126) num3,trms3,chi3
          write(15,126) num3,trms3,chi3
          if(itomo.ne.2) write(38,126) num3,trms3,chi3
126       format('for rays turning beneath the interfce:'/
     +           'RMS traveltime residual for ',i8,' rays: ',f12.5/
     +         '                   Normalized chi-squared: ',f12.4/)
        end if
      end if
c
      if(ioop.eq.1) then
        write(6,6175) 100.*soop/noop,noop,100.*msoop/noop
6175    format(/'Average out-of-plane distance (%):',f7.2/
     +          'Number of raypaths used:          ',i7/
     +          'Average max. out-of-plane distance (%):',f7.2/)
        do i=1,100
           if(oophist(i).gt.0) then
             l1oop(i)=alog10(1.*oophist(i))
           else
             l1oop(i)=0.
           end if
           if(moophist(i).gt.0) then
              l2oop(i)=alog10(1.*moophist(i))
           else
             l2oop(i)=0.
           end if
        enddo
        write(56,6185) (i-1,oophist(i),l1oop(i),i=1,100)
        write(57,6185) (i-1,moophist(i),l2oop(i),i=1,100)
6185    format(2i10,f10.3)
      end if
c
      if(itomo.eq.5) then
c
        if(nlc.gt.0) then 
          write(41) (dkernel(j),j=1,nwrite) 
          write(42) (row(j),j=1,nwrite)
          write(43) (column(j),j=1,nwrite)
        endif 
c
        write(io,3135) nl
        write(44,*) nl,num2
3135    format('number of non-zero elements in data kernel:',i10/)
      end if
c
      if(iplots.eq.1) then
        call empty
        call plotnd(1)
      end if
c
      if(itomo.eq.1) then 
        rewind(17)
        rewind(18)
        rewind(29)
        do 3131 k=1,inz 
           do 3131 j=1,iny 
              write(17) (pert1(i,j,k),i=1,inx) 
              write(18) (pert2(i,j,k),i=1,inx) 
3131          write(29) (numray(i,j,k),i=1,inx)
      end if
c
      if(itomo.eq.5) then 
        rewind(29)
        do 3132 k=1,inz 
           do 3132 j=1,iny 
3132          write(29) (numray(i,j,k),i=1,inx)
      end if
c
      if(itomo.eq.3) then 
        rnorm=0.
        if(i2d.ne.1) then
          do 4110 k=2,inz-1
             do 4110 j=2,iny-1
                do 4110 i=2,inx-1
                   pert2(i,j,k)=(pert1(i,j,k)+pert1(i,j,k-1)+
     +                       pert1(i,j-1,k)+pert1(i,j-1,k-1)+
     +                       pert1(i-1,j,k)+pert1(i-1,j,k-1)+
     +                  pert1(i-1,j-1,k)+pert1(i-1,j-1,k-1))/8.
                   if(pert2(i,j,k).gt.rnorm) rnorm=pert2(i,j,k)
4110      continue
        else
          do 4130 k=2,inz
             do 4130 i=2,inx
                pert2(i,1,k)=0.
                do 4140 j=1,iny
4140               pert2(i,1,k)=pert1(i,j,k)+pert1(i,j,k-1)+
     +                       pert1(i-1,j,k)+pert1(i-1,j,k-1)
                pert2(i,1,k)=pert2(i,1,k)/(4.*float(iny))
                do 4150 j=1,iny
4150               pert2(i,j,k)=pert2(i,1,k)
                if(pert2(i,1,k).gt.rnorm) rnorm=pert2(i,1,k)
4130      continue
        end if

        if(rnorm.eq.0.) rnorm=1.
c
        if(i2d.ne.1) then
          do 4120 k=2,inz
             do 4120 j=2,iny
                do 4120 i=2,inx
4120               pert1(i,j,k)=pert2(i,j,k)/rnorm
        else
          do 4160 k=2,inz
             do 4160 j=1,iny
                do 4160 i=2,inx
4160               pert1(i,j,k)=pert2(i,j,k)/rnorm
        end if
c
        rewind(17)
        do 3150 k=1,inz 
           do 3150 j=1,iny 
3150          write(17) (pert1(i,j,k),i=1,inx) 
      end if
c
      stop        
      end         
