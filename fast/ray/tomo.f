c
c     version 1.0  Nov 1993
c
c     tomography related routines for RAY
c                 
c        patched for g77 compilation
c
c        Drew Brenders (brenders@geoladm.geol.queensu.ca)
c
c     ----------------------------------------------------------------
c                 
      subroutine backproj(tobs,uobs,tcalc,interfce,iflag)
c                 
c     calculate the slowness perturbation for each model cell sampled 
c     by the current ray path using back projection
c                 
      include 'ray.par'
      include 'ray.com'
c
      real length
      integer iback(nray),ci(nray,3),iipt(nray)
c
      length=0.
      iflag=0
      nipts=0
      ixl=0
      iyl=0
      izl=0
      nback=0
c
c     resample raypath onto inverse grid
c
      do 50 i=1,npts-1
         xmid=(xray(i)+xray(i+1))/2.
         ymid=(yray(i)+yray(i+1))/2.
         zmid=(zray(i)+zray(i+1))/2.
c
         ix=min(int((xmid-xmin)/xii+1),inx)
         if(ny.ne.1) then
           iy=min(int((ymid-ymin)/yii+1),iny)
         else
           iy=1
         end if
         iz=min(int((zmid-zmin)/zii+1),inz)
c
         if(ix.ne.ixl.or.iy.ne.iyl.or.iz.ne.izl) then
           nipts=nipts+1
           lseg(nipts)=sqrt((xray(i)-xray(i+1))**2+
     +            (yray(i)-yray(i+1))**2+(zray(i)-zray(i+1))**2)
           ixl=ix
           iyl=iy
           izl=iz
           ci(nipts,1)=ix
           ci(nipts,2)=iy
           ci(nipts,3)=iz
         else
           lseg(nipts)=lseg(nipts)+sqrt((xray(i)-xray(i+1))**2+
     +           (yray(i)-yray(i+1))**2+(zray(i)-zray(i+1))**2)
         end if
         iipt(i)=nipts
50    continue
c
c     determine effective length of raypath
c
      if(interfce.eq.0) then
        do 10 i=1,nipts
           length=length+lseg(i)
10      continue
      else
        do 31 i=1,nipts
31         iback(i)=1
c
        do 30 i=1,npts-1
           z1=bndinterp(xray(i),yray(i),cell(i,1),cell(i,2),1)
           if(z1.gt.zray(i)) iback(iipt(i))=0
           z2=bndinterp(xray(i+1),yray(i+1),cell(i,1),cell(i,2),1)
           if(z2.gt.zray(i+1)) iback(iipt(i))=0
30      continue
c
        do 70 i=1,nipts
           if(iback(i).ne.0) length=length+lseg(i)
70      continue
      end if
c
      if(length.eq.0.) then
        iflag=1
        return
      end if
c
      deltat=tobs-tcalc
c
      if(imethod.eq.0) then
        sperturbation=deltat/length
      else
        sum=0.
        if(interfce.eq.0) then
          do 80 i=1,nipts
             sum=sum+lseg(i)*imodel(ci(i,1),ci(i,2),ci(i,3))
80        continue
        else
          do 90 i=1,nipts
             if(iback(i).eq.1) 
     +         sum=sum+lseg(i)*imodel(ci(i,1),ci(i,2),ci(i,3))
90        continue
        end if
        const=deltat/sum
      end if
c
      if(iunc.eq.0) then
        unc=1.
      else
        unc=uobs
      end if
c
      if(interfce.eq.0) then
        do 20 i=1,nipts
           if(imethod.eq.1) sperturbation=
     +       const*imodel(ci(i,1),ci(i,2),ci(i,3))
           pert1(ci(i,1),ci(i,2),ci(i,3))=
     +     pert1(ci(i,1),ci(i,2),ci(i,3))+sperturbation/unc
           pert2(ci(i,1),ci(i,2),ci(i,3))=
     +     pert2(ci(i,1),ci(i,2),ci(i,3))+1./unc
           numray(ci(i,1),ci(i,2),ci(i,3))=
     +     numray(ci(i,1),ci(i,2),ci(i,3))+one
20      continue
      else
        do 40 i=1,nipts
           if(iback(i).eq.1) then
             if(imethod.eq.1) sperturbation=
     +         const*imodel(ci(i,1),ci(i,2),ci(i,3))
             pert1(ci(i,1),ci(i,2),ci(i,3))=
     +       pert1(ci(i,1),ci(i,2),ci(i,3))+sperturbation/unc
             pert2(ci(i,1),ci(i,2),ci(i,3))=
     +       pert2(ci(i,1),ci(i,2),ci(i,3))+1./unc
             numray(ci(i,1),ci(i,2),ci(i,3))=
     +       numray(ci(i,1),ci(i,2),ci(i,3))+one
           end if
40      continue
      end if
c
      return
      end
c                 
c     ----------------------------------------------------------------
c                 
      subroutine d_time(tobs,uobs,tcalc,interfce,iflag)
c                 
c     calculate the time perturbation associated with a known
c     slowness perturbation 
c                 
      include 'ray.par'
      include 'ray.com'
c
      real length
      integer iback(nray),ci(nray,3),iipt(nray)
c
      length=0.
      iflag=0
      nipts=0
      ixl=0
      iyl=0
      izl=0
      nback=0
c
c     resample raypath onto inverse grid
c
      do 50 i=1,npts-1
         xmid=(xray(i)+xray(i+1))/2.
         ymid=(yray(i)+yray(i+1))/2.
         zmid=(zray(i)+zray(i+1))/2.
c
         ix=min(int((xmid-xmin)/xii+1),inx)
         if(ny.ne.1) then
           iy=min(int((ymid-ymin)/yii+1),iny)
         else
           iy=1
         end if
         iz=min(int((zmid-zmin)/zii+1),inz)
c
         if(ix.ne.ixl.or.iy.ne.iyl.or.iz.ne.izl) then
           nipts=nipts+1
           lseg(nipts)=sqrt((xray(i)-xray(i+1))**2+
     +            (yray(i)-yray(i+1))**2+(zray(i)-zray(i+1))**2)
           ixl=ix
           iyl=iy
           izl=iz
           ci(nipts,1)=ix
           ci(nipts,2)=iy
           ci(nipts,3)=iz
         else
           lseg(nipts)=lseg(nipts)+sqrt((xray(i)-xray(i+1))**2+
     +           (yray(i)-yray(i+1))**2+(zray(i)-zray(i+1))**2)
         end if
         iipt(i)=nipts
50    continue
c
c     determine effective length of raypath
c
      if(interfce.eq.0) then
        do 10 i=1,nipts
           length=length+lseg(i)
10      continue
      else
        do 31 i=1,nipts
31         iback(i)=1
c
        do 30 i=1,npts-1
           z1=bndinterp(xray(i),yray(i),cell(i,1),cell(i,2),1)
           if(z1.gt.zray(i)) iback(iipt(i))=0
           z2=bndinterp(xray(i+1),yray(i+1),cell(i,1),cell(i,2),1)
           if(z2.gt.zray(i+1)) iback(iipt(i))=0
30      continue
c
        do 70 i=1,nipts
           if(iback(i).ne.0) length=length+lseg(i)
70      continue
      end if
c
      if(length.eq.0.) then
        iflag=1
        return
      end if
c
      deltat=tobs-tcalc
      if(iunc.eq.0) then
        unc=1.
      else
        unc=uobs
      end if
      sum=0.
c
      if(interfce.eq.0) then
        do 20 i=1,nipts
           sum=sum+pert1(ci(i,1),ci(i,2),ci(i,3))*lseg(i)
20      continue
      else
        do 40 i=1,nipts
           if(iback(i).eq.1) then
             sum=sum+pert1(ci(i,1),ci(i,2),ci(i,3))*lseg(i)
           end if
40      continue
      end if
c
      write(28) deltat/uobs,sum/uobs
c
      return
      end
c                 
c     ----------------------------------------------------------------
c                 
      subroutine kernel(tobs,uobs,uave,tcalc,interfce,iflag)
c                 
c     calculate the raypath length in each model cell of inverse
c     grid for current ray path
c                 
      include 'ray.par'
      include 'ray.com'
c
      real length(nximax*nyimax*nzimax)
c
      iflag=1
c
      do 20 i=1,nm
20       length(i)=0.
c
      do 10 i=1,npts-1
c
         if(interfce.eq.1) then
           z1=bndinterp(xray(i),yray(i),cell(i,1),cell(i,2),1)
           if(z1.gt.zray(i)) go to 10
           z2=bndinterp(xray(i+1),yray(i+1),cell(i,1),cell(i,2),1)
           if(z2.gt.zray(i+1)) go to 10
         end if
         iflag=0
c
         xl=((xray(i)-xray(i+1))**2+(yray(i)-yray(i+1))**2+
     +       (zray(i)-zray(i+1))**2)**.5
c
         xmid=(xray(i)+xray(i+1))/2.
         ymid=(yray(i)+yray(i+1))/2.
         zmid=(zray(i)+zray(i+1))/2.
c
         ix=min(int((xmid-xmin)/xii+1),inx)
         if(ny.ne.1) then
           iy=min(int((ymid-ymin)/yii+1),iny)
         else
           iy=1
         end if
         iz=min(int((zmid-zmin)/zii+1),inz)
c
         im=ix+(iy-1)*inx+(iz-1)*inx*iny
         length(im)=length(im)+xl
         numray(ix,iy,iz)=1
10    continue
c
      if(iflag.eq.1) return
c
      nk=nk+1
      uw=uave/uobs
c
      nlc=0
      do 40 i=1,nm
         if(length(i).gt.0.) then
           nlc=nlc+1
c           dkernel(nlc)=uw*length(i)
c           row(nlc)=nk
c           column(nlc)=i
c           print*,'kernel',nlc,column(nlc)
c           if(nlc.eq.nwrite) then
c             print*,'WRITE!'
c             write(41) (dkernel(j),j=1,nwrite)
cc             write(42) (row(j),j=1,nwrite)
c             write(43) (column(j),j=1,nwrite)
c             nlc=0
c           endif
c           print*,'Kernel',nl,uw*length(i),nk,i
           write(41) uw*length(i)
           write(42) nk
           write(43) i
             nlc=0



           nl=nl+1
         end if
40    continue
c
      if(istype.ne.2) write(40) uw*(tobs-tcalc)
c
      return
      end
c        
c     ----------------------------------------------------------------
c        
      subroutine resolution(ix0,iy0,iz0,ix1,iy1,iz1)
c      
c     calculate the resolution kernel (point spread function) centred
c     at the node (ix0,iy0,iz0)
c      
      include 'ray.par'
      include 'ray.com'
c
      real length
      length=0.
      wlength=0.
c
c     determine effective length of raypath
c
      do 10 i=1,npts-1
         lseg(i)=((xray(i)-xray(i+1))**2+(yray(i)-yray(i+1))**2+
     +           (zray(i)-zray(i+1))**2)**.5
         length=length+lseg(i)
         if((cell(i,1).eq.ix0.or.cell(i,1).eq.ix1).and.
     +      (cell(i,2).eq.iy0.or.cell(i,2).eq.iy1).and.
     +      (cell(i,3).eq.iz0.or.cell(i,3).eq.iz1)) 
     +      wlength=wlength+lseg(i)
10    continue
c
      if(length.eq.0.) return
      res=wlength/length
c      
      do 20 i=1,npts-1
         pert1(cell(i,1),cell(i,2),cell(i,3))=
     +   pert1(cell(i,1),cell(i,2),cell(i,3))+res
20    continue
c
      return
      end
