c
      open(unit=11, file='rec.out', form='unformatted')
      nrec=0
c
      write(6,15)
15    format('enter shot position (x,y,z)')
      read(5,*) xs,ys,zs
      write(11) xs,ys,zs,0.,0.,-1
c
100   write(6,1)
1     format(/
     +'enter 1 for linear array, 2 for 2D array (0 to stop)')
      read(5,*) iarray
      if(iarray.le.0) go to 999
c
      if(iarray.ne.1) then
      write(6,2)
2     format(/
     +'enter xmin,xmax,ymin,ymax,zmin,zmax for receiver locations')
      read(5,*) xrmin,xrmax,yrmin,yrmax,zrmin,zrmax
      write(6,3)
3     format('enter receiver spacing in x,y,z directions')
      read(5,*) xinc,yinc,zinc
c
      nrx=(xrmax-xrmin)/xinc+1
      nry=(yrmax-yrmin)/yinc+1
      nrz=(zrmax-zrmin)/zinc+1
c
      do 10 i=1,nrx
         xr=xrmin+float(i-1)*xinc
         do 10 j=1,nry
            yr=yrmin+float(j-1)*yinc
            do 10 k=1,nrz
               zr=zrmin+float(k-1)*zinc
               write(11) xr,yr,zr,0.,.05,1
10             nrec=nrec+1
      else
      write(6,4) 
4     format(/'enter point (x,y,z) for point 1')
      read(5,*) xrmin,yrmin,zrmin
      write(6,5) 
5     format(/'enter point (x,y,z) for point 2')
      read(5,*) xrmax,yrmax,zrmax
      write(6,6)
6     format('enter receiver spacing')
      read(5,*) rinc
c
      dist=((xrmin-xrmax)**2+(yrmin-yrmax)**2+(zrmin-zrmax)**2)**.5
      nr=nint(dist/rinc)+1
      xinc=(xrmax-xrmin)/float(nr-1)
      yinc=(yrmax-yrmin)/float(nr-1)
      zinc=(zrmax-zrmin)/float(nr-1)
c
      do 20 i=1,nr
         xr=xrmin+float(i-1)*xinc
         yr=yrmin+float(i-1)*yinc
         zr=zrmin+float(i-1)*zinc
         write(11) xr,yr,zr,0.,.05,1
20       nrec=nrec+1
      end if
c
      go to 100
c
999   write(6,25) nrec
25    format(/'number of receivers: ',i10/)
c
      stop
      end
