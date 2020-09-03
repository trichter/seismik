c        patched for g77 compilation
c
c        Drew Brenders (brenders@geoladm.geol.queensu.ca)
c
c
      character*72 ifname,ofname
c
      write(*, fmt="(/'Enter input file name')")
      read(5,85) ifname
85    format(a72)
      write(*, fmt="(/'Enter output file name')")
      read(5,85) ofname
c
      open(unit=11, file=ifname, form='unformatted', status='old')
      open(unit=12, file=ofname)
      n=0
c
100   read(11, end=999) x,y,z,t,u,i
      write(12,5) x,y,z,t,u,i
5     format(5f10.3,i3)
      n=n+1
      go to 100
c
999   write(6,1) n
1     format(/'number of lines in file: ',i10/)
c
      stop
      end
