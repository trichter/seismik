c
      character*72 ifname,ofname
      integer*2 v(601)
c
      write(*, fmt="(/'Enter input file name')")
      read(5,85) ifname
85    format(a72)
      write(*, fmt="(/'Enter output file name')")
      read(5,85) ofname
c
      open(unit=11, file=ifname, form='unformatted', status='old')
      open(unit=12, file=ofname )
      n=0
c     
      do 101 j=1,25
100      read(11) (v(i), i=1,601) 
         do 102 i=1,601
            write(12,*) v(i)
102      continue
         write(*,*) v(601)

         n=n+1
c      go to 100
101   end do






999   write(6,1) n
1     format(/'number of lines in file: ',i10/)
c
      stop
      end
