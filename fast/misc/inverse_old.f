c
c     version 1.0  Oct 1998
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |               ********  I N V E R S E  *****                 |   
c     |                                                              |
c     |   Calculate the slowness update using LSQR regularization    |
c     |                      below a boundary                        |
c     |                                                              |
c     |                 Written by Colin A. Zelt                     |   
c     |                                                              |
c     |                Dept of Geology & Geophysics                  |   
c     |                      Rice University                         |
c     |                  Houston, TX  77005-1892                     |   
c     |                                                              |
c     ----------------------------------------------------------------
c
c     smallest and smoothest (or flattest) perturbation with edge 
c     constraint and with variable weighting of perturbation and 
c     horizontal smoothness (flatness) with depth
c
c     I/O units:
c
c        10 -- input:  inverse 3D data grid
c
c        11 -- output:  input parameters
c
c        35 -- input:  parameters describing size of inverse 3D grid
c
c
c     ----------------------------------------------------------------
c
c
      include 'ray.par'
c
c     ndmax = number of data
c     nnzmax = number of non-zero elements in matrix L from ray
c
c    ***********************************************************
c    * If you change nnzmax and/or ndmax, make sure you change *
c    * it in the subroutines avpu and atupv below!!!!!!!!!!!!! *
c    ***********************************************************
c
      parameter(nnzmax=200000, ndmax=3400)
      parameter(nsconst=4*nximax*nyimax*nzimax,
     +          nszero=9*nximax*nyimax*nzimax,
     +          nmmax=nximax*nyimax*nzimax)
c
      integer row(nnzmax+nszero),column(nnzmax+nszero),
     +        ireg(nximax,nyimax),florsm(nzimax)/nzimax*2/
      real lambda,imodel(nximax,nyimax,nzimax),dm(nmmax),
     +     dum1(nmmax),dum2(nmmax),data(ndmax+nsconst),
     +     anz(nnzmax+nszero),ibmodel(nximax,nyimax,nzimax),
     +     smwz(nzimax)/nzimax*-1./,spwz(nzimax)/nzimax*-1./
c
      common /blk1/ anz,row,column,nnz
c
      namelist /invpar/ sz,error,itmax,iupdate,kstart,alpha,
     +                  interface,smwz,spwz,sedge,florsm
c
      sedge=10.
      interface=0
      alpha=0.95
      iupdate=1
      error=.000000001
      sz=0.2
      itmax=50
      kstart=0
c
      write(6,335)
335   format('INVERSE: ',
     +'calculate slowness update using LSQR and regularization')
      open(50, file='nowrite', status='old', err=9999, iostat=ioflag)
9999  if(ioflag.gt.0) then
        io=0
      else
        io=51
      end if
c
      open(48, file='stop.in', status='old')
      read(48,*) iistop
      if(iistop.gt.0) stop
c
      open(11, file='i.in',status='old')
c
      read(11,invpar)
c
c     check the alpha value
c
      if(alpha.lt.0.or.alhpa.gt.1) then
        write(0,11)
11      format(/'***  0 <= alpha <= 1 required  ***')
        stop
      end if
      sm=alpha
      sp=1.-alpha
      if(alpha.lt..00001) then
        sm=0. 
        sp=1.
      end if
      if(alpha.gt..99999) then
        sm=1. 
        sp=0.
      end if
      if(sz.lt..00001) sz=0.
c
      open(13, file='delta.time', form='unformatted', status='old')
      open(12, file='data.kernel', form='unformatted', status='old')
      open(14, file='row.kernel', form='unformatted', status='old')
      open(15, file='column.kernel', form='unformatted', status='old')
      open(44, file='nzero.kernel', status='old')
      open(19, file='sl.pert', form='unformatted')
      open(10, file='slow.mod', form='unformatted', status='old')
      open(16, file='ref_slow.mod', form='unformatted', status='old')
      open(35, file='inv.header', status='old')
      open(45, file='lambda', status='old')
c
      read(45,*) lambda
c
      read(35,*) nx,ny,nz
c
      if(nx.gt.nximax.or.ny.gt.nyimax.or.nz.gt.nzimax) then
        write(0,9)
9       format(/'***  model is too large  ***'/)
        stop
      end if
c
      nmodel=nx*ny*nz
      write(io,25) nx,ny,nz,nmodel
25    format(/'model dimensions:'/
     +        '-----------------'/
     +        'number of model cells in each direction (x,y,z):',3i4/
     +        'total number of model cells:                    ',i12)
c
      if(interface.eq.1) then
        open(22, file='ireg.cell', form='unformatted',status='old')
        do j=1,ny
           read(22) (ireg(i,j),i=1,nx)
        enddo
      end if
c
      do 10 k=1,nz
         do 10 j=1,ny
         read(10) (imodel(i,j,k),i=1,nx)
10       read(16) (ibmodel(i,j,k),i=1,nx)
c
      read(44,*) nk,nobs
c
      if(nk.gt.nnzmax) then
        write(0,*) '***  increase nnzmax  ***'
        stop
      end if
c
      nr=nk/nwrite
      nl=nk-nr*nwrite
      if(nl.gt.0) nr=nr+1
c
      do 20 j=1,nr
         is=(j-1)*nwrite
         read(12) (anz(is+i),i=1,nwrite)
         read(14) (row(is+i),i=1,nwrite)
         read(15) (column(is+i),i=1,nwrite)
20    continue
c
      do 30 i=1,nobs
30       read(13) data(i)
c
      sx=lambda*sm
      szz=sx*sz
      sx4=4.*sx
      sz2=2.*szz
      spert=lambda*sp
      sedge=sedge*lambda
      nconst=nobs
      nnz=nk
c
      if(sm.gt.0.) then
c
      do 400 k=1,nz
         if(smwz(k).ge.0.) then
           smw=smwz(k)
         else
           smw=1.
         end if
         if(smw.eq.0.) go to 400
         if(florsm(k).ne.1) then
c
c          horizontal smoothness
c
           do 40 j=2,ny-1
              do 40 i=2,nx-1
           if(interface.eq.1) then
           imax=max(ireg(i,j),ireg(i-1,j),ireg(i+1,j),
     +              ireg(i,j-1),ireg(i,j+1))
           else
           imax=0
           end if
           if(imax+kstart.le.k.or.interface.eq.0) then    
           im1=i+(j-1)*nx+(k-1)*nx*ny
           im2=im1+1
           im3=im1-1
           im4=i+(j-2)*nx+(k-1)*nx*ny
           im5=i+(j)*nx+(k-1)*nx*ny
           nconst=nconst+1
c
           sn=imodel(i,j,k)
           nnz=nnz+1
           anz(nnz)=smw*sx4/sn
           row(nnz)=nconst
           column(nnz)=im1
           nnz=nnz+1
           anz(nnz)=-smw*sx/sn
           row(nnz)=nconst
           column(nnz)=im2
           nnz=nnz+1
           anz(nnz)=-smw*sx/sn
           row(nnz)=nconst
           column(nnz)=im3
           nnz=nnz+1
           anz(nnz)=-smw*sx/sn
           row(nnz)=nconst
           column(nnz)=im4
           nnz=nnz+1
           anz(nnz)=-smw*sx/sn
           row(nnz)=nconst
           column(nnz)=im5
c
           data(nconst)=
     +       smw*(-sx4*(imodel(i,j,k)-ibmodel(i,j,k))+
     +       sx*(imodel(i+1,j,k)-ibmodel(i+1,j,k))+
     +       sx*(imodel(i-1,j,k)-ibmodel(i-1,j,k))+
     +       sx*(imodel(i,j+1,k)-ibmodel(i,j+1,k))+
     +       sx*(imodel(i,j-1,k)-ibmodel(i,j-1,k)))/sn
           end if
40         continue
         else
c
c          x-direction horizontal flatness
c
           do 41 j=1,ny
              do 41 i=1,nx-1
           if(interface.eq.1) then
             imax=max(ireg(i,j),ireg(i+1,j))
           else
             imax=0
           end if
           if(imax+kstart.le.k.or.interface.eq.0) then
           im1=i+(j-1)*nx+(k-1)*nx*ny
           im2=im1+1
           nconst=nconst+1
c
           sn=(imodel(i,j,k)+imodel(i+1,j,k))/2.
           nnz=nnz+1
           anz(nnz)=smw*sx/sn
           row(nnz)=nconst
           column(nnz)=im1
           nnz=nnz+1
           anz(nnz)=-smw*sx/sn
           row(nnz)=nconst
           column(nnz)=im2
c
           data(nconst)=
     +       smw*sx*(-(imodel(i,j,k)-ibmodel(i,j,k))+
     +       (imodel(i+1,j,k)-ibmodel(i+1,j,k)))/sn
           end if
41         continue
c
c          y-direction horizontal flatness
c
           do 42 j=1,ny-1
              do 42 i=1,nx
           if(interface.eq.1) then
             imax=max(ireg(i,j),ireg(i,j+1))
           else
             imax=0
           end if
           if(imax+kstart.le.k.or.interface.eq.0) then
           im1=i+(j-1)*nx+(k-1)*nx*ny
           im2=i+(j)*nx+(k-1)*nx*ny
           nconst=nconst+1
c
           sn=(imodel(i,j,k)+imodel(i,j+1,k))/2.
           nnz=nnz+1
           anz(nnz)=smw*sx/sn
           row(nnz)=nconst
           column(nnz)=im1
           nnz=nnz+1
           anz(nnz)=-smw*sx/sn
           row(nnz)=nconst
           column(nnz)=im2
c
           data(nconst)=
     +       smw*sx*(-(imodel(i,j,k)-ibmodel(i,j,k))+
     +       (imodel(i,j+1,k)-ibmodel(i,j+1,k)))/sn
           end if
42         continue
         end if
400   continue
c
c     vertical smoothness 
c
      if(sz.gt.0.) then
c
      do 60 k=2,nz-1
         do 60 j=1,ny
            do 60 i=1,nx
         if(interface.eq.1) then
           imax=ireg(i,j)
         else
           imax=0
         end if
         if(imax+kstart.le.k-1.or.interface.eq.0) then
         im1=i+(j-1)*nx+(k-2)*nx*ny
         im2=i+(j-1)*nx+(k-1)*nx*ny
         im3=i+(j-1)*nx+(k)*nx*ny
         nconst=nconst+1
c
         sn=imodel(i,j,k)
         nnz=nnz+1
         anz(nnz)=szz/sn
         row(nnz)=nconst
         column(nnz)=im1
         nnz=nnz+1
         anz(nnz)=-sz2/sn
         row(nnz)=nconst
         column(nnz)=im2
         nnz=nnz+1
         anz(nnz)=szz/sn
         row(nnz)=nconst
         column(nnz)=im3
c
         data(nconst)=(-szz*(imodel(i,j,k-1)-ibmodel(i,j,k-1))+
     +                 sz2*(imodel(i,j,k)-ibmodel(i,j,k)) 
     +                 -szz*(imodel(i,j,k+1)-ibmodel(i,j,k+1)))/sn
         end if
60    continue
c
      end if
c
      end if
c
c     model perturbation
c
      do 1100 k=1,nz
         if(spwz(k).ge.0.) then
           spw=spwz(k)
         else
           spw=1.
         end if
         do 110 j=1,ny
            do 110 i=1,nx
               if(interface.eq.1) then
                 imax=ireg(i,j)
               else
                 imax=0
               end if
               if(imax+kstart.le.k.or.interface.eq.0) then
                 im1=i+(j-1)*nx+(k-1)*nx*ny
c
                 if(i.eq.1.or.j.eq.1.or.k.eq.1.or.
     +              i.eq.nx.or.j.eq.ny.or.k.eq.nz) then
                   if(sedge.gt.0.) then
                     nconst=nconst+1
                     sn=imodel(i,j,k)
                     nnz=nnz+1
                     anz(nnz)=sedge/sn
                     row(nnz)=nconst
                     column(nnz)=im1
                     data(nconst)=-anz(nnz)*(imodel(i,j,k)-
     +                                   ibmodel(i,j,k))
                   end if
                 else
                   if(sp.gt.0.) then
                     nconst=nconst+1
                     sn=imodel(i,j,k)
                     nnz=nnz+1
                     anz(nnz)=spw*spert/sn
                     row(nnz)=nconst
                     column(nnz)=im1
                     data(nconst)=-anz(nnz)*(imodel(i,j,k)-
     +                                   ibmodel(i,j,k))
                   end if
                 end if
               end if
110      continue
1100  continue
c
      if(nnz.gt.nnzmax+nszero) then
        write(0,*) '***  increase nnzmax or nsconst  ***'
        stop
      end if
      if(nconst.gt.ndmax+nsconst) then
        write(0,*) '***  increase ndmax or nsconst  ***'
        stop
      end if
c
      write(io,*) 
      write(io,45) nmodel,nconst,nnz
45    format(/'number of model parameters:     ',i10/
     +        'number of constraint equations: ',i10/
     +        'number of non-zero elements:    ',i10)
      write(io,35) float(nnz)/(float(nconst)*float(nmodel))*100.
35    format('% non-zero elements: ',f10.6)
c
c     solve the system using LSQR
c
      call lsqr(nconst,nmodel,dm,data,dum1,dum2,itmax,it,
     +          error,iupdate,io,r1,rlast)
c
c     output slowness update
c
      do 190 k=1,nz
         do 190 j=1,ny
            do 210 i=1,nx
            im=i+(j-1)*nx+(k-1)*nx*ny
210         imodel(i,j,k)=dm(im)
            write(19) (imodel(i,j,k),i=1,nx)
190   continue
c
      open(38, file='log.file',status='old',form='formatted')
155   read(38,55,end=98) alog
55    format(a1)
      go to 155
c
98    backspace(38)
      write(38,45) nmodel,nconst,nnz
      write(38,35) float(nnz)/(float(nconst)*float(nmodel))*100.
      write(38,65) it,sz,alpha,sedge/lambda,kstart
65    format('number of LSQR iterations: ',i10/
     +'sz:',e10.4,'  alpha:',e10.4,'  sedge:',e10.4,'  kstart:',i6)
      write(38,75) r1,rlast
75    format('first and last relative error:',2f10.6)
c
      stop
      end
c
c     ----------------------------------------------------------------
c
      subroutine lsqr(m,n,x,u,v,w,itmax,it,error,iupdate,io,r1,rlast)
c
      real x(n),u(m),v(n),w(n)
c
      if(iupdate.eq.1) then
        open(29, file='error')
        write(29,*) error
        close(29)
      end if
c
      do i=1,n
         x(i)=0.
         v(i)=0.
      end do
c
      call normlz(m,u,beta)
c
      b1=beta
c
      call atupv(m,n,u,v)
      call normlz(n,v,alfa)
c
      rhobar=alfa
      phibar=beta
c
      do i=1,n
         w(i)=v(i)
      end do
c
      do iter=1,itmax
         a=-alfa
         do i=1,m
            u(i)=a*u(i)
         end do
c
         call avpu(m,n,u,v)
         call normlz(m,u,beta)
c
         b=-beta
c
         do i=1,n
            v(i)=b*v(i)
         end do
c
         call atupv(m,n,u,v)
         call normlz(n,v,alfa)
c
         rho=sqrt(rhobar*rhobar+beta*beta)
         c=rhobar/rho
         s=beta/rho
         teta=s*alfa
         rhobar=-c*alfa
         phi=c*phibar
         phibar=s*phibar
         t1=phi/rho
         t2=-teta/rho
c
         do i=1,n
            x(i)=t1*w(i)+x(i)
            w(i)=t2*w(i)+v(i)
         end do
c
         r=phibar/b1
c
         if(iter.eq.1)  then
           write(io,5) iter,phibar,r,0.
           r1=r
         else
           write(io,5) iter,phibar,r,rprev-r
           rlast=r
         end if
5        format(i5,e15.6,2f10.6)
c
         if(iupdate.eq.1) then
           open(29, file='error', status='old')
           read(29,*) error
           close(29)
         end if
c
         if(abs(r-rprev).lt.error) then
           write(io,*) 'number of iterations in LSQR: ',iter
           it=iter
           return
         end if
c
         rprev=r
      end do
c
      it=itmax
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine normlz(n,x,s)
c
      real x(n)
c
      s=0.
c
      do i=1,n
         s=s+x(i)*x(i)
      end do
c
      s=sqrt(s)
      ss=1./s
c
      do i=1,n
         x(i)=x(i)*ss
      end do
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine avpu(m,n,u,v)
c
      include 'ray.par'
c
      parameter(nnzmax=200000, ndmax=3400)
      parameter(nsconst=4*nximax*nyimax*nzimax,
     +          nszero=9*nximax*nyimax*nzimax,
     +          nmmax=nximax*nyimax*nzimax)
c
      integer row(nnzmax+nszero),column(nnzmax+nszero)
      real anz(nnzmax+nszero),u(m),v(n),sum(ndmax+nsconst)
c
      common /blk1/ anz,row,column,nnzero
c
      do i=1,m
         sum(i)=0.
      end do
c
      do i=1,nnzero
         sum(row(i))=sum(row(i))+anz(i)*v(column(i))
      end do
c
      do i=1,m
         u(i)=u(i)+sum(i)
      end do
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine atupv(m,n,u,v)
c
      include 'ray.par'
c
      parameter(nnzmax=200000, ndmax=3400)
      parameter(nsconst=4*nximax*nyimax*nzimax,
     +          nszero=9*nximax*nyimax*nzimax,
     +          nmmax=nximax*nyimax*nzimax)
c
      integer row(nnzmax+nszero),column(nnzmax+nszero)
      real anz(nnzmax+nszero),u(m),v(n),sum(nmmax)
c
      common /blk1/ anz,row,column,nnzero
c
      do j=1,n
         sum(j)=0.
      end do
c
      do i=1,nnzero
         sum(column(i))=sum(column(i))+anz(i)*u(row(i))
      end do
c
      do j=1,n
         v(j)=v(j)+sum(j)
      end do
c
      return
      end
