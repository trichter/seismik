c
c     version 1.1  Apr 1995
c
c     model routine for FD
c
c     ----------------------------------------------------------------
c
      subroutine model
c
c     read in and setup the velocity model
c
      include 'fd.par'
      include 'fd.com'
c
      integer*2 veli(nxmax)
c
c     read in model in 3D binary format
c
      open(34, file='for.header', status='old')
c
c     read in model dimensions
c
      read(34,115) xmin,xmax,ymin,ymax,zmin,zmax,size,nx,ny,nz
115   format(7f10.3,3i10)
c
      if(nx.gt.nxmax.or.nz.gt.nzmax.or.ny.gt.nymax) then
        write(0,1)
1       format(/'***  model is too large  ***'/)
        stop
      end if
c
      open(23, file='vel.mod', form='unformatted',
     +     status='old')
c
c     convert the velocity array to velocity/size
c
      vmult=.001/size
      do 50 k=1,nz
         do 50 j=1,ny
           read(23) (veli(i),i=1,nx)
           do 51 i=1,nx
51            vel(i,j,k)=veli(i)*vmult
50    continue
c
      sidelimit(1)=1
      sidelimit(2)=nx
      sidelimit(3)=1
      sidelimit(4)=ny
      sidelimit(5)=1
      sidelimit(6)=nz
      nodestotal=nx*ny*nz
      size2=size**2
c
      open(50, file='nowrite', status='old', err=9999, iostat=ioflag)
9999  if(ioflag.gt.0) then
      write(io,25) xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,nodestotal
25    format(/'model dimensions:'/
     +        '-----------------'/
     +        'xmin, xmax: ',2f10.3/'ymin, ymax: ',2f10.3/,
     +        'zmin, zmax: ',2f10.3/
     +        'number of model nodes in each direction (x,y,z):',3i5/
     +        'total number of model nodes:                    ',i15)
      write(io,35) size
35    format('node spacing: ',f7.3,' km')
      end if
c
      if(i2d.eq.1) then
        write(io,55)
55      format(/'>>>  running in 2D mode  <<<')
        nodestotal=nx*nz
      end if
c
      return
      end         
c
c     ----------------------------------------------------------------
c
      subroutine find(xp,yp,zp,ic,jc,kc)
c
c     find which model cell the point (xp,yp,zp) is in
c
      include 'fd.par'
      include 'fd.com'
c
      ic=int((xp-xmin)/size)+1
      jc=int((yp-ymin)/size)+1
      kc=int((zp-zmin)/size)+1
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine findnode(xs,ys,zs,ixs,iys,izs)
c
c     find closest model node to the source
c
      include 'fd.par'
      include 'fd.com'
c
      ixs=nint((xs-xmin)/size)+1
      iys=nint((ys-ymin)/size)+1
      izs=nint((zs-zmin)/size)+1
c
      return
      end
c
c     ----------------------------------------------------------------
c
      function velocity(xp,yp,zp,ixc,iyc,izc)
c
c     calculate the velocity at the point (xp,yp,zp)
c
      include 'fd.par'
      include 'fd.com'
c
      v1=vel(ixc,iyc,izc)
      v2=vel(ixc+1,iyc,izc)

      v5=vel(ixc,iyc,izc+1)
      v6=vel(ixc+1,iyc,izc+1)

c
      xfrac=xp-(float(ixc-1)*size+xmin)
      yfrac=yp-(float(iyc-1)*size+ymin)
      zfrac=zp-(float(izc-1)*size+zmin)
c
      if(ny.eq.1) then
        v3=v1
        v4=v2
        v7=v5
        v8=v6
        yfrac=0.
      else
        v3=vel(ixc,iyc+1,izc)
        v4=vel(ixc+1,iyc+1,izc)
        v7=vel(ixc,iyc+1,izc+1)
        v8=vel(ixc+1,iyc+1,izc+1)
      end if
c
      vleft= (size*(v3-v1)*yfrac+size*(v5-v1)*zfrac+
     +       (v1-v3+v7-v5)*yfrac*zfrac)/size2+v1 
      vright=(size*(v4-v2)*yfrac+size*(v6-v2)*zfrac+
     +       (v2-v4+v8-v6)*yfrac*zfrac)/size2+v2 
c
      velocity=(size-xfrac)*vleft+xfrac*vright
c
      return
      end
c
c     ----------------------------------------------------------------
c
      function vaverage(ix1,iy1,iz1,ix2,iy2,iz2)
c
c     calculate the average velocity/size inside a cell with opposite 
c     corners at points 1 and 2 using a simple average of the 8
c     values at the corners
c
      include 'fd.par'
      include 'fd.com'
c
      vaverage=(vel(ix1,iy1,iz1)+vel(ix1,iy1,iz2)+
     +          vel(ix1,iy2,iz1)+vel(ix1,iy2,iz2)+
     +          vel(ix2,iy1,iz1)+vel(ix2,iy1,iz2)+
     +          vel(ix2,iy2,iz1)+vel(ix2,iy2,iz2))*.125
c
      return
      end
