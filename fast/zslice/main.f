c
c     version 1.1  Apr 1995
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |            **********    Z S L I C E  *********              |
c     |                                                              |
c     |             Plot 2D slices of a 3D data volume               |
c     |                                                              |
c     |                 Written by Colin A. Zelt                     |
c     |                                                              |
c     |                   Bullard Laboratories                       |
c     |                  University of Cambridge                     |
c     |                  Cambridge, UK  CB3 0EZ                      |
c     |                                                              |
c     ----------------------------------------------------------------
c
c
c        patched for g77 compilation
c
c        Drew Brenders (brenders@geoladm.geol.queensu.ca)
c 
c     I/O units:  
c                 
c        10 -- input:  program input parameters
c                 
c        11 -- input:  input parameters for panning
c                 
c        12 -- input:  file names and related parameters
c                 
c        19 -- output: all Calcomp plot calls
c
c        21 -- input: interface
c
c        35 -- input:  parameters describing size of 3D grid
c
c        36 -- input:  3D data
c                 
c
c     ----------------------------------------------------------------
c                 
c 
      program main
c
      include 'zslice.par'
      include 'zslice.com'                 
c 
      real xplot(nx2dmax),yplot(ny2dmax),xpl(2),ypl(2),
     +     inter(nxmax,nymax),intplot(nx2dmax)
      character plane*1,xlab*6,ylab*6,hfile*25,dfile*25,kclass*1,ph*1,
     +          label*80,ifile*25,reply*1
      integer*2 timei(nxmax),inti(nxmax)
c
      namelist /pltpar/ iplot,iroute,hfile,
     +                  iseg,xwndow,ywndow,ibcol,ifcol,colour,iscale,
     +                  symht,igrid,anht,ndecia,angap,icolor,itype,
     +                  xsc,int2,ndecic,icolneg,icolpos,ipan,
     +                  interface,intcol,ifile,irect,ncntr,
     +                  power,ictable,dc
c
      namelist /axepar/ iaxlab,albht,xorig,yorig,
     +                  xmm,ymm,zmm
c
c     initialize parameters
c
      dc=0.
      ncntr=21
      intcol=-1
      interface=0
      itype=0
      ipan=0
      icolneg=1
      icolpos=3
      undef=999.999
      iscan=0
      icolor=0
      ndecic=1
      int2=0
      igrid=0
      symht=1.
      isep=0
      xorig=15.
      yorig=-1.
      iroute=1
      iaxlab=1
c
      open(10, file='s.in', status='old')
      open(12, file='sfile.in', status='old')
c      
c     read in program control parameters from unit 10
c                 
      read(10,pltpar)
      read(10,axepar)
c
c     open I/O units
c
      if(iplot.eq.0) iplot=-1
      if(iplot.eq.2) iplot=0
      if(iplot.le.0) open(19, file='p.out')
      open(35, file=hfile, status='old')
      if(ipan.eq.1) then
        open(11, file='zslice.pan', status='old')
        iin=11
      else
        iin=5
      end if
c
c     read in data
c
      read(35,115) xmin,xmax,ymin,ymax,zmin,zmax,size,nx,ny,nz
115   format(7f10.3,3i10)
c
      nodestotal=nx*ny*nz
      write(6,25) xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,nodestotal
25    format(/'model dimensions:'/
     +        '-----------------'/
     +        'xmin, xmax: ',2f10.3/'ymin, ymax: ',2f10.3/,
     +        'zmin, zmax: ',2f10.3/
     +        'number of model nodes in each direction (x,y,z):',3i4/
     +        'total number of model nodes:                    ',i12)
      write(6,35) size
35    format('node spacing: ',f7.3,' km')
c
      if(nx.gt.nxmax.or.ny.gt.nymax.or.nz.gt.nzmax) then
        write(6,3)
3       format(/'***  model is too large  ***'/)
        stop
      end if
c
9999  rewind(12)
      i=1
      read(12,165)
165   format(' ')
333   read(12,*,end=222) dfile
      write(6,125) i,dfile
125   format('   ',i3,' - ',a25)
      if(mod(i,2).eq.0) write(6,165)
      i=i+1
      go to 333
222   rewind(12)
      read(12,165)
146   write(6,135)
135   format(/'enter data file number (0 to stop)')
      read(5,*,err=146) idata
      if(idata.eq.0) then
        if(iplots.eq.1) call plotnd(ipan)
        stop
      end if
      do 147 i=1,abs(idata)
147      read(12,*) dfile,int2,cmin,cmax,cinc,scalef
c
      nk=0
      write(6,145) dfile
145   format('data file: ',a25)
c
      if(itype.ne.1) then
c
        open(36, file=dfile, form='unformatted', 
     +       status='old')
c
        if(int2.ne.1) then
          do 129 k=1,nz
             do 130 j=1,ny
130             read(36,end=1399) (time(i,j,k),i=1,nx)
             nk=nk+1
129       continue
          if(scalef.ne.1.) then
            do 131 k=1,nz
               do 131 j=1,ny
                  do 131 i=1,nx
131                    time(i,j,k)=time(i,j,k)*scalef
          end if
        else
          do 128 k=1,nz
             do 132 j=1,ny
                read(36,end=1399) (timei(i),i=1,nx)
                do 133 i=1,nx
                   time(i,j,k)=float(timei(i))*scalef
133                if(timei(i).eq.-32768) time(i,j,k)=undef
132          continue
             nk=nk+1
128       continue
        end if
      else
c
        if(int2.ne.1) then
c
          open(36, file=dfile, form='unformatted',
     +       access='direct', recl=4*nx, status='unknown')
c
          do 136 k=1,nz
             do 137 j=1,ny
c137             read(36,rec=((k-1)*ny+j),end=1399) 
137             read(36,end=1399) (time(i,j,k),i=1,nx)
             nk=nk+1
136       continue
          if(scalef.ne.1.) then
            do 138 k=1,nz
               do 138 j=1,ny
                  do 138 i=1,nx
138                  time(i,j,k)=time(i,j,k)*scalef
          end if
        else
c
          open(36, file=dfile, form='unformatted',
     +       access='direct', recl=2*nx, status='unknown')
c
          do 139 k=1,nz
             do 141 j=1,ny
c                read(36,rec=((k-1)*ny+j),end=1399) 
                read(36,end=1399) (timei(i),i=1,nx)
                do 142 i=1,nx
                   time(i,j,k)=float(timei(i))*scalef
142                if(timei(i).eq.-32768) time(i,j,k)=undef
141          continue
             nk=nk+1
139       continue
        end if
      end if
c
      go to 1400
c
1399  if(nk.eq.1) then
        write(6,155) nk
155     format(/'>>> data file contains ',i3,' surface'/)
      else
        write(6,156) nk
156     format(/'>>> data file contains ',i3,' surfaces'/)
      end if
      do 134 k=nk+1,nz
         do 134 j=1,ny
            do 134 i=1,nx
134            time(i,j,k)=time(i,j,nk)
c
1400  close(36)
      if(interface.eq.1) then
        open(21, file=ifile, form='unformatted', status='old')
c        do 5120 j=1,ny
c5120       read(21) (inti(i),i=1,nx)
        do 5140 k=1,nz
        do 5140 j=1,ny
c           read(21) (inti(i),i=1,nx)
            read(21,end=1499) (inter(i,j),i=1,nx)
           do 5150 i=1,nx
5150          inter(i,j)=inter(i,j)/1000.
c5150          inter(i,j)=inti(i)/1000.
5140    continue
1499    if(intcol.lt.0) intcol=ifcol
        close(21)
      end if
c
c     calculate scale of each plot axis
c
      xscalea=(xmax-xmin)/xmm
      yscalea=(ymax-ymin)/ymm
      zscalea=-(zmax-zmin)/zmm
      if(yorig.lt.0.) yorig=xorig
      xinc=(.999999*xmax-xmin)/float(ninterp-1)
      yinc=(.999999*ymax-ymin)/float(ninterp-1)
      zinc=(.999999*zmax-zmin)/float(ninterp-1)
      size2=size**2
c
      tamin=1.e20       
      tamax=-1.e20
      ixmin=1
      iymin=1
      izmin=1
      ixmax=1
      iymax=1
      izmax=1
      do 150 i=1,nx
         do 150 j=1,ny
            do 150 k=1,nz
               if(time(i,j,k).lt.tamin) then
                 tamin=time(i,j,k)
                 ixmin=i
                 iymin=j
                 izmin=k
               end if
               if(time(i,j,k).gt.tamax.and.
     +         abs(time(i,j,k)-undef).gt..001) then
                 tamax=time(i,j,k)
                 ixmax=i
                 iymax=j
                 izmax=k
               end if
150   continue
c
      ximin=xmin+float(ixmin-1)*size
      yimin=ymin+float(iymin-1)*size
      zimin=zmin+float(izmin-1)*size
      ximax=xmin+float(ixmax-1)*size
      yimax=ymin+float(iymax-1)*size
      zimax=zmin+float(izmax-1)*size
      write(6,2) ximin,yimin,zimin,tamin,ximax,yimax,zimax,tamax
2     format('min data value at (',f8.3,',',f8.3,',',f8.3,'):',f12.5/
     +       'max data value at (',f8.3,',',f8.3,',',f8.3,'):',f12.5)
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
        colour(6)=7
        colour(7)=8
        colour(8)=9
        colour(9)=10
        colour(10)=11
        colour(11)=12
        colour(12)=13
        colour(13)=14
        colour(14)=15
        ncol=14
      end if
c
      if(cmin.gt.999999.) cmin=tamin
      if(cmax.lt.-999999.) cmax=tamax
      if(cinc.le.0.) then
        ncntrs=21
        cinc=(cmax-cmin)/20.
      else
        ncntrs=nint((cmax-cmin)/cinc)+1
      end if
      cinc2=cinc/2.
      iflag1=0
c                 
111   continue
c
      if(iscan.eq.1) then
        nscan=nscan+1
        position=pmin+float(nscan-1)*pinc
        if(position.gt.pmax) iscan=0
        go to 1010
      end if
      iiflag=0
c
795   write(6,1)
1     format(/
     +'enter normal direction to slice plane (x,y or z),',
     +' position (km),'/
     +'and plot type (1,2,3) (enter position=9999 for scanning'
     +,' or s,0,0 to stop)')
      read(iin,*,err=795,end=9999) plane,position,imeth
      if(imeth.lt.0) then
        iout=1
      else
        iout=0
      end if
      imeth=iabs(imeth)
      if(plane.eq.'s') go to 9999
      if(plane.eq.'v'.and.iflag1.eq.1) go to 2221
      if(plane.ne.'x'.and.plane.ne.'y'.and.plane.ne.'z')
     +   go to 795
      if(imeth.lt.1.or.imeth.gt.3) go to 795
c
      if(position.eq.9999.) then
        write(6,715)
715     format('enter spacing of slice planes for scanning')
        read(iin,*) pinc
        iscan=1
        nscan=0
        if(plane.eq.'x') then
          pmin=xmin
          pmax=xmax
        end if
        if(plane.eq.'y') then
          pmin=ymin
          pmax=ymax
        end if
        if(plane.eq.'z') then
          pmin=zmin
          pmax=zmax
        end if
        go to 111
      end if
c
2221  if(plane.eq.'v') then
c
2222    call get_event(xc,yc,keyval,kclass)
c
        if(kclass.eq.'v') then
1111      write(6,505)
505       format(
     +    'enter x,y,z coordinate for data value (0,0,0 to stop)')
          read(iin,*,err=1111) xp,yp,zp
          ixp=int((xp-xmin)/size)+1
          iyp=int((yp-ymin)/size)+1
          izp=int((zp-zmin)/size)+1
          write(6,515) timeinterp(xp,yp,zp,ixp,iyp,izp)
515       format('data value: ',f12.5)
          if(xp.eq.0.and.yp.eq.0.and.zp.eq.0.) go to 111
          go to 1111
        end if
c
        if(xc.ge.xorig.and.xc.le.xorig+xlen.and.
     +     yc.ge.yorig.and.yc.le.yorig+ylen) then
          if(ph.eq.'x') then
            xp=xpos
            yp=(xc-xorig)*xscale+xpmin
            zp=(yc-yorig)*yscale+ypmax
          end if
          if(ph.eq.'y') then
            yp=ypos
            xp=(xc-xorig)*xscale+xpmin
            zp=(yc-yorig)*yscale+ypmax
          end if
          if(ph.eq.'z') then
            zp=zpos
            xp=(xc-xorig)*xscale+xpmin
            yp=(yc-yorig)*yscale+ypmax
          end if
          ixp=int((xp-xmin)/size)+1
          iyp=int((yp-ymin)/size)+1
          izp=int((zp-zmin)/size)+1
          write(6,525) xp,yp,zp,timeinterp(xp,yp,zp,ixp,iyp,izp)
525       format('(x,y,z): ',3f10.3,'    data value: ',f12.5)
          go to 2222
        else
          go to 111
        end if
c
      end if
c
1010  continue
c
      if(plane.eq.'z') then
        if(position.lt.zmin.or.position.gt.zmax) go to 111
        ph=plane
        kpos=nint((position-zmin)/size)+1
        npx=nx
        npy=ny
        xpmin=xmin
        xpmax=xmax
        ypmin=ymin
        ypmax=ymax
        xlen=xmm
        ylen=ymm
        xscale=xscalea
        yscale=-yscalea
        yadj=ymax
        xlab='X (km)'
        ylab='Y (km)'
        zpos=position
        tamin=1.e20           
        tamax=-1.e20
        ixmin=1
        iymin=1
        ixmax=1
        iymax=1
        iflagm=0
c       trms=0.
        do 10 i=1,nx
           do 10 j=1,ny
              if(time(i,j,kpos).lt.tamin) then
                 tamin=time(i,j,kpos)
                 ixmin=i
                 iymin=j
               end if
               if(time(i,j,kpos).gt.tamax.and.
     +         abs(time(i,j,kpos)-undef).gt..001) then
                 tamax=time(i,j,kpos)
                 ixmax=i
                 iymax=j
                 iflagm=1
               end if
               if(iout.eq.1) write(37,*) xmin+float(i-1)*size,
     +                       ymax-float(j-1)*size,time(i,j,kpos)
c              trms=trms+time(i,j,kpos)**2
10             tplot(i,j)=time(i,j,kpos)
c       trms=sqrt(trms/float(nx*ny))
        ximin=xmin+float(ixmin-1)*size
        yimin=ymin+float(iymin-1)*size
        ximax=xmin+float(ixmax-1)*size
        yimax=ymin+float(iymax-1)*size
        if(iflagm.eq.1) write(6,4) ximin,yimin,tamin,ximax,yimax,tamax
4       format('min data value of slice at (',f8.3,',',f8.3,'):',f12.5/
     +         'max data value of slice at (',f8.3,',',f8.3,'):',f12.5)
c       write(0,445) trms
445     format('trms: ',f10.5)
      end if
c
      if(plane.eq.'y') then
        if(position.lt.ymin.or.position.gt.ymax) go to 111
        ph=plane
        jpos=nint((position-ymin)/size)+1
        npx=nx
        npy=nz
        xpmin=xmin
        xpmax=xmax
        ypmin=zmin
        ypmax=zmax
        xlen=xmm
        ylen=zmm
        xscale=xscalea
        yscale=zscalea
        yadj=zmax
        xlab='X (km)'
        ylab='Z (km)'
        ypos=position
        tamin=1.e20           
        tamax=-1.e20
        ixmin=1
        izmin=1
        ixmax=1
        izmax=1
        iflagm=0
c       trms=0.
        do 20 i=1,nx
           do 20 k=1,nz
              if(time(i,jpos,k).lt.tamin) then
                 tamin=time(i,jpos,k)
                 ixmin=i
                 izmin=k
               end if
               if(time(i,jpos,k).gt.tamax.and.
     +         abs(time(i,jpos,k)-undef).gt..001) then
                 tamax=time(i,jpos,k)
                 ixmax=i
                 izmax=k
                 iflagm=1
               end if
               if(iout.eq.1) write(37,*) xmin+float(i-1)*size,
     +                       -zmin-float(k-1)*size,time(i,jpos,k)
c              trms=trms+time(i,jpos,k)**2
20             tplot(i,k)=time(i,jpos,k)
c       trms=sqrt(trms/float(nx*nz))
        ximin=xmin+float(ixmin-1)*size
        zimin=zmin+float(izmin-1)*size
        ximax=xmin+float(ixmax-1)*size
        zimax=zmin+float(izmax-1)*size
        if(iflagm.eq.1) write(6,4) ximin,zimin,tamin,ximax,zimax,tamax
        if(interface.eq.1) then
          do 210 i=1,npx
210          intplot(i)=inter(i,jpos)
          iiflag=1
        end if
c       write(0,445) trms
      end if
c
      if(plane.eq.'x') then
        if(position.lt.xmin.or.position.gt.xmax) go to 111
        ph=plane
        ipos=nint((position-xmin)/size)+1
        npx=ny
        npy=nz
        xpmin=ymin
        xpmax=ymax
        ypmin=zmin
        ypmax=zmax
        xlen=ymm
        ylen=zmm
        xscale=yscalea
        yscale=zscalea
        yadj=zmax
        xlab='Y (km)'
        ylab='Z (km)'
        xpos=position
        tamin=1.e20           
        tamax=-1.e20
        iymin=1
        izmin=1
        iymax=1
        izmax=1
        iflagm=0
        do 30 j=1,ny
           do 30 k=1,nz
              if(time(ipos,j,k).lt.tamin) then
                 tamin=time(ipos,j,k)
                 iymin=j
                 izmin=k
               end if
               if(time(ipos,j,k).gt.tamax.and.
     +         abs(time(ipos,j,k)-undef).gt..001) then
                 tamax=time(ipos,j,k)
                 iymax=j
                 izmax=k
                 iflagm=1
               end if
               if(iout.eq.1) write(37,*) ymin+float(j-1)*size,
     +                       -zmin-float(k-1)*size,time(ipos,j,k)
30             tplot(j,k)=time(ipos,j,k)
        yimin=ymin+float(iymin-1)*size
        zimin=zmin+float(izmin-1)*size
        yimax=ymin+float(iymax-1)*size
        zimax=zmin+float(izmax-1)*size
        if(iflagm.eq.1) write(6,4) yimin,zimin,tamin,yimax,zimax,tamax
        if(interface.eq.1) then
          do 310 j=1,npx
310          intplot(j)=inter(ipos,j)
          iiflag=1
        end if
      end if
c
      sizex=size/xscale
      sizey=size/yscale
      xtmin=-999999.
      ytmin=-999999.
      xtmax=-999999.
      ytmax=-999999.
      ndecix=-2
      ndeciy=-2
      ntickx=-1
      nticky=-1
c
      if(iplots.eq.0) then
        call plots(xwndow,ywndow,iroute)
        call segmnt(1)
        iplots=1
      else
        if(ipan.ne.1) then
          call erase
          call segmnt(1)
        end if
      end if
c
      if(ipan.eq.1) then
66      write(6,65)
65      format(/'enter xorig,yorig,ixlab,iylab,label')
        read(iin,*,end=9999, err=66) xorig,yorig,ixlab,iylab,label
        ixlab=2*ixlab-1
        iylab=2*iylab-1
        xxlabel=xlabel+xorig
        yylabel=ylabel+yorig
      else
        ixlab=1
        iylab=1
      end if
c
      if(iaxlab.eq.1) then
        call pcolor(ifcol)
        call axtick(xpmin,xpmax,xtmin,xtmax,ntickx,ndecix)
        call axis(xorig,yorig,xpmin,xpmax,xlen,xscale,0.,1,
     +       xtmin,xtmax,ntickx,ndecix,xlab,6,ixlab*albht)
        call axtick(ypmin,ypmax,ytmin,ytmax,nticky,ndeciy)
        call axis(xorig,yorig,ypmin,ypmax,ylen,yscale,90.,1,
     +       ytmin,ytmax,nticky,ndeciy,ylab,6,iylab*albht)
      end if
c
      call pcolor(ifcol)
      call box(xorig,yorig,xorig+xlen,yorig+ylen)
c
      do 80 i=1,npx 
         xplot(i)=xpmin+float(i-1)*size
80       xplot(i)=(xplot(i)-xpmin)/xscale+xorig
      do 90 j=1,npy 
         yplot(j)=ypmin+float(j-1)*size
90       yplot(j)=(yplot(j)-yadj)/yscale+yorig
c
      if(igrid.eq.1) then
        call pcolor(2)
        do 60 i=2,npx-1
           call plot(xplot(i),yorig,3)
60         call plot(xplot(i),yorig+ylen,2)
        do 70 j=2,npy-1
           call plot(xorig,yplot(j),3)
70         call plot(xorig+xlen,yplot(j),2)
        call pcolor(ifcol)
      end if
c
      if(imeth.ge.2) then
           do 500 i=1,npx-1,ndecic
              do 600 j=1,npy-1,ndecic
                 if(abs(max(tplot(i,j),tplot(i+1,j),tplot(i,j+1),
     +                 tplot(i+1,j+1))-undef).lt.001) go to 600
                 tave=(tplot(i,j)+tplot(i+1,j)+tplot(i,j+1)+
     +                 tplot(i+1,j+1))/4.
                 if(icolor.eq.0) then
                   if(tave.lt.cmin) go to 600
                   if(tave.gt.cmax) then
                      dotht=size/xscale
                   else
                      dotht=(tave-cmin)/(cmax-cmin)*size/xscale
                   end if
                   icol=ifcol
                 else
                   if(tave.gt.-cinc2.and.tave.lt.cinc2) go to 600
                   if(tave.lt.dc) then
                     if(tave.lt.cmin) then
                        dotht=size/xscale
                     else
                        dotht=(abs(tave/cmin))*size/xscale
                     end if
                     icol=colour(icolneg)
                   else
                     if(tave.gt.cmax) then
                        dotht=size/xscale
                     else
                        dotht=(tave/cmax)*size/xscale
                     end if
                     icol=colour(icolpos)
                   end if
                 end if
                 call dot(xplot(i)+.5*size/xscale-dotht/2.,
     +                yplot(j)+.5*size/yscale+dotht/2.,dotht,icol)
600           continue
500        continue
      end if
c
      if(imeth.eq.0) then
        do 100 ic=1,ncntrs
           icol=mod(ic,ncol)
           tcntr=cmin+float(ic-1)*cinc
           do 100 i=1,npx,ndecic
             do 100 j=1,npy-1,ndecic
               if((tplot(i,j).lt.tcntr.and.tplot(i,j+1).ge.tcntr).or.
     +            (tplot(i,j).ge.tcntr.and.tplot(i,j+1).lt.tcntr)) then
                  ypp=sizey/(tplot(i,j+1)-tplot(i,j))*
     +                (tcntr-tplot(i,j))+yplot(j) 
                  call dot(xplot(i),ypp,symht,colour(icol))
               end if
100     continue
c
        do 200 ic=1,ncntrs
           icol=mod(ic,ncol)
           tcntr=cmin+float(ic-1)*cinc
           do 200 j=1,npy,ndecic
             do 200 i=1,npx-1,ndecic
               if((tplot(i,j).lt.tcntr.and.tplot(i+1,j).ge.tcntr).or.
     +            (tplot(i,j).ge.tcntr.and.tplot(i+1,j).lt.tcntr)) then
                  xpp=sizex/(tplot(i+1,j)-tplot(i,j))*
     +                (tcntr-tplot(i,j))+xplot(i)
                  call dot(xpp,yplot(j),symht,colour(icol))
               end if
200     continue
      end if
c
      if(imeth.eq.1.or.imeth.eq.3) then
        nn2=0
        ifirst=0
        if(idata.lt.0) then
          cmin=tamin
          cmax=tamax
          ncntrs=ncntr
          cinc=(cmax-cmin)/float(ncntrs-1)
          cinc2=cinc/2.
        end if
        do 300 ic=1,ncntrs
           tc=cmin+float(ic-1)*cinc
           if(tc.lt.tamin.or.tc.gt.tamax) go to 300
           if(icolor.eq.0) then
             call pcolor(ifcol)
           else
             if(tc.lt.dc) then
               call pcolor(icolneg+1)
             else
               call pcolor(icolpos+1)
             end if
           end if
           do 400 i=1,npx-1,ndecic
              do 400 j=1,npy-1,ndecic
                 t1=tplot(i,j)
                 t2=tplot(i+1,j)
                 t3=tplot(i,j+1)
                 t4=tplot(i+1,j+1)
                 np=0
                 if((t1.lt.tc.and.t2.ge.tc).or.(t1.ge.tc.and.t2.lt.tc)) 
     +           then
                   np=np+1
                   if(np.gt.2) then
                     nn2=nn2+1
                     go to 400
                   end if
                   xpl(np)=sizex/(t2-t1)*(tc-t1)+xplot(i)
                   ypl(np)=yplot(j)
                 end if
                 if((t3.lt.tc.and.t4.ge.tc).or.(t3.ge.tc.and.t4.lt.tc)) 
     +           then
                   np=np+1
                   if(np.gt.2) then
                     nn2=nn2+1
                     go to 400
                   end if
                   xpl(np)=sizex/(t4-t3)*(tc-t3)+xplot(i)
                   ypl(np)=yplot(j+1)
                 end if
                 if((t1.lt.tc.and.t3.ge.tc).or.(t1.ge.tc.and.t3.lt.tc)) 
     +           then
                   np=np+1
                   if(np.gt.2) then
                     nn2=nn2+1
                     go to 400
                   end if
                   xpl(np)=xplot(i)
                   ypl(np)=sizey/(t3-t1)*(tc-t1)+yplot(j)
                 end if
                 if((t2.lt.tc.and.t4.ge.tc).or.(t2.ge.tc.and.t4.lt.tc)) 
     +           then
                   np=np+1
                   if(np.gt.2) then
                     nn2=nn2+1
                     go to 400
                   end if
                   xpl(np)=xplot(i+1)
                   ypl(np)=sizey/(t4-t2)*(tc-t2)+yplot(j)
                 end if
                 if(np.lt.2) then
                   if(nn2.eq.1) nn2=nn2+1
                   go to 400
                 end if
                 if(ifirst.eq.0) then
                   tfirst=tc
                   ifirst=1
                 end if
                 tlast=tc
                 call plot(xpl(1),ypl(1),3)
                 call plot(xpl(2),ypl(2),2)
400        continue
300     continue
        if(ifirst.eq.1) write(6,75) tfirst,tlast
75      format('first and last contour values: ',2f10.3)
        if(nn2.ne.0) write(6,55) nn2
55      format(
     +  'number of cells with 1 or more than 2 intersection points: ',
     +  i10)
      end if
c
      if(iiflag.eq.1) then
        do 910 i=1,npx
910        intplot(i)=(intplot(i)-yadj)/yscale+yorig
        call pcolor(intcol)
        call line(xplot,intplot,npx)
      end if
c
      call empty
c
      if(iscan.eq.1) then
        write(6,705) position
705     format('position: ',f10.3,' - enter <CR> for next slice,'
     +,         ' s to stop scan')
        read(iin,101) reply
c        read(*,101) reply
101     format(a1)
        if(reply.eq.'s') iscan=0
      end if
c
      iflag1=1
      go to 111
c
      end         
