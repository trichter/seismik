c        patched for g77 compilation
c
c        Drew Brenders (brenders@geoladm.geol.queensu.ca)
c
c
      character*72 if1name,if2name
c
      write(*, fmt="(/'Enter input file name')")
      read(5,85) if1name
      write(*, fmt="(/'Enter input file name')")
      read(5,85) if2name
85    format(a72)
c
      open(unit=11, file=if1name, form='unformatted', status='old')
      open(unit=12, file=if2name, form='unformatted', status='old')
      open(unit=13, file='rec_diff.out')
      n=0
      dmax=0.
      sum=0.
      sum1=0.
      sum2=0.
c
100   read(11, end=999) x,y,z,t1,u,i
      read(12, end=999) x,y,z,t2,u,i
      if(i.eq.-1) then
        write(13,5) x,y,z,t1,u,i
5       format(5f10.3,i3)
      else
        write(13,5) x,y,z,t2-t1,u,i
        sum=sum+t2-t1
        sum1=sum1+abs(t2-t1)
        sum2=sum2+(t2-t1)**2
        if(abs(t2-t1).gt.dmax) then
          dmax=abs(t2-t1)
          xmax=x
          ymax=y
          zmax=z
        end if
        n=n+1
      end if
      go to 100
c
999   write(6,1) n
1     format(/'number of times in file: ',i10/)
c
      ave=sum/float(n)
      ave1=sum1/float(n)
      ave2=sqrt(sum2/float(n))
c
      write(6,2) ave,ave1,ave2
2     format(/'averge time difference:           ',f10.4/
     +        'average absolute time difference: ',f10.4/
     +        'RMS time difference:              ',f10.4)
      write(6,3) dmax,xmax,ymax,zmax
3     format(/'maximum time difference: ',f10.4/
     +        'at (x,y,z):              ',3f10.3/)
c
      stop
      end
