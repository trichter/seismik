c        patched for g77 compilation
c
c        Drew Brenders (brenders@geoladm.geol.queensu.ca)
c
c
      character*72 ifname,ofname
      num=0
c
      write(6,5)
5     format(/
     +'Enter seed (positive integer)')
      read(5,*) iseed
      rnum=gasdev(-iseed)
c
1000  write(*, fmt="(/'Enter input file name')")
      read(5,85) ifname
      if(ifname.eq.'') then
        write(6,2) num
2       format(/'total number of lines: ',i12)
        stop
      end if
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
      if(i.ge.0) then
        rnum=u*gasdev(iseed)
        t=t+rnum
      end if
c
      write(12) x,y,z,t,u,i
      n=n+1
      go to 100
c
999   write(6,1) n
1     format(/'number of lines in file: ',i10/)
c
      num=num+n
      close(11)
      close(12)
c
      go to 1000
c
      end
c
c     -----------------------------------------------------------------
c
      function gasdev(idum)
c
c     returns a normally distributed number with zero mean and unit
c     variance
c
      data iset/0/
      if(iset.eq.0) then
1       v1=2.*ran1(idum)-1.
        v2=2.*ran1(idum)-1.
        r=v1**2+v2**2
        if(r.ge.1.) go to 1
        fac=sqrt(-2.*log(r)/r)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      end if
      return
      end
c
c     -----------------------------------------------------------------
c
      function ran1(idum)
c
c     returns a uniform random number between 0.0 and 1.0. Set IDUM to
c     any negative value to initialize or reinitialize the sequence.
c
      dimension r(97)
      parameter (m1=259200, ia1=7141, ic1=54773, rm1=1./m1)
      parameter (m2=134456, ia2=8121, ic2=28411, rm2=1./m2)
      parameter (m3=243000, ia3=4561, ic3=51349)
      data iff/0/
      if(idum.lt.0.or.iff.eq.0) then
        iff=1
        ix1=mod(ic1-idum,m1)
        ix1=mod(ia1*ix1+ic1,m1)
        ix2=mod(ix1,m2)
        ix1=mod(ia1*ix1+ic1,m1)
        ix3=mod(ix1,m3)
        do 11 j=1,97
           ix1=mod(ia1*ix1+ic1,m1)
           ix2=mod(ia2*ix2+ic2,m2)
           r(j)=(float(ix1)+float(ix2)*rm2)*rm1
11      continue
        idum=1
      end if
      ix1=mod(ia1*ix1+ic1,m1)
      ix2=mod(ia2*ix2+ic2,m2)
      ix3=mod(ia3*ix3+ic3,m3)
      j=1+(97*ix3)/m3
      if(j.gt.97.or.j.lt.1) pause
      ran1=r(j)
      r(j)=(float(ix1)+float(ix2)*rm2)*rm1
      return
      end
