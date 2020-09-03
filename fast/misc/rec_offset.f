c        patched for g77 compilation
c
c        Drew Brenders (brenders@geoladm.geol.queensu.ca)
c

      character*72 ifname,ofname
c
      write(6,5)
5     format(/
     +'Enter minimum and maximum offset (km)')
      read(5,*) xmin,xmax
c
1000  write(*, fmt="(/'Enter input file name')")
      read(5,85) ifname
      if(ifname.eq.'') stop
c
85    format(a72)
      write(*, fmt="(/'Enter output file name (default is input)')")
      read(5,85) ofname
      if(ofname.eq.'') ofname=ifname
c
      open(unit=11, file=ifname, form='unformatted', status='old')
      open(unit=12, file=ofname, form='unformatted')
c
      n=0
c
100   read(11,end=999) x,y,z,t,u,i
c
      if(i.eq.-1) then
        xs=x
        ys=y
        zs=z
        write(12) x,y,z,t,u,i
        n=n+1
      else
        offset=sqrt((x-xs)**2+(y-ys)**2+(z-zs)**2)
        if(offset.ge.xmin.and.offset.le.xmax) then
          write(12) x,y,z,t,u,i
          n=n+1
        end if
      end if
c
      go to 100
c
999   write(6,1) n
1     format(/'number of lines in file: ',i10/)
c
      close(11)
      close(12)
c
      go to 1000
c
      end
