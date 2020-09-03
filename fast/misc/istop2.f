c
c     version 1.0  Jul 1996
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |              **********  I S T O P  **********               |
c     |                                                              |
c     |  Check chi**2 values to see if iterations should be stopped  |
c     |                                                              |
c     |                 Written by Colin A. Zelt                     |   
c     |                                                              |
c     |                Dept of Geology & Geophysics                  |   
c     |                      Rice University                         |
c     |                  Houston, TX  77005-1892                     |   
c     |                                                              |
c     ----------------------------------------------------------------
c
c
c     I/O units:
c
c        11 -- input: chi.values
c
c        13 -- input: stop.in
c
c
c     ----------------------------------------------------------------
c
      include 'ray.par'
c
      character*1 alog
      real lambda
      tol=.005
c
      write(6,335)
335   format('ISTOP: check to see if iterations should be stopped')
      open(50, file='nowrite', status='old', err=9999, iostat=ioflag)
9999  if(ioflag.gt.0) then
        io=0
      else
        io=51
      end if
c
      open(11, file='chi.values', status='old')
      open(13, file='stop.in', status='old')
      open(14, file='lambda', status='old')
c
      read(14,*) lambda
c
c     read(11,1) 
c     read(11,1) 
      read(11,*) 
      read(11,*) 
      read(11,*) 
1     format(' ')
c
      chimin=1.e20
c
      i=1
100   read(11,5,end=999) chi
5     format(43x,f11.5)

      if((chi.gt.chimin.and.chimin.gt.1.).or.abs(chi-1.).lt.tol) then
        write(13,*) 1
c
        open(38, file='log.file',status='old',form='formatted')
155     read(38,55,end=98) alog
55      format(a1)
        go to 155
c
98      write(6,45) 
        backspace(38)
        write(38,45) 
45      format(/'***  stop iterating  ***')
c
        stop
      end if
c
      if(chi.lt.chimin) chimin=chi 
      i=i+1
      go to 100
c
999   stop
      end
