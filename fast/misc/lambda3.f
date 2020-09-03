c
c     version 1.0  Apr 1996
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |           **********  L A M B D A 2  **********              |
c     |                                                              |
c     |           Calculate the trade-off parameter lambda           |
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
c        11 -- output:  lambda value
c
c        13 -- input:  input parameters        
c
c     ----------------------------------------------------------------
c
c
      parameter(nmax=100)
c
      real lambda0,lambda,plamb(nmax),pchi(nmax),jumfac,lambl1,lambg1
      integer chione
c
      namelist /lampar/ redfac,lambda0,jumfac
c
      jumfac=3.
      lambda0=100.
      redfac=2.
c
      write(6,335)
335   format('LAMBDA2: calculate trade-off parameter')
      open(50, file='nowrite', status='old', err=9999, iostat=ioflag)
9999  if(ioflag.gt.0) then
        io=0
      else
        io=51
      end if
c
      open(49, file='stop.in', status='old')
      read(49,*) iistop
      if(iistop.gt.0) stop
c
      open(14, file='lambda', status='old')
      open(13, file='l.in', status='old')
c
      read(13,lampar) 
c
      read(14,*) lambda
c
      if(lambda.eq.0.) then
        lambda=lambda0
      else
        if(lambda.lt.0.) then
          lambda=abs(lambda)
        else
c
          open(11, file='chi.values', status='old')
          open(12, file='lambda.values', status='old')
c
c         read(11,1)
c         read(11,1)
1         format(' ')
           read(11,*)
           read(11,*)

          i=1
100       read(11,5,end=99) pchi(i)
5         format(43x,f11.5)
          read(12,15,end=99) plamb(i)
15        format(17x,e15.6)
          i=i+1
          go to 100

99        np=i-1
c
          lambg1=0.
          lambl1=0.
          if(np.gt.1) then
            chig1=1.e20
            chil1=0.
            do 10 i=1,np
               if(pchi(i).gt.1..and.pchi(i).lt.chig1) then
                 chig1=pchi(i)
                 lambg1=plamb(i)
               end if
               if(pchi(i).lt.1..and.pchi(i).gt.chil1) then
                 chil1=pchi(i)
                 lambl1=plamb(i)
               end if

10          continue
          end if
c
          if(lambl1.gt.0..and.np.gt.1) then
            if(lambg1.eq.0.) then
              chig1=0.
              do 20 i=1,np
                 if(pchi(i).gt.chig1.and.plamb(i).ne.lambl1) then
                   chig1=pchi(i)
                   lambg1=plamb(i)
                 end if
20            continue
            end if
            slope=(lambg1-lambl1)/(chig1-chil1)
            b=lambg1-slope*chig1
            lambda=slope+b  
            lambda=min(lambda,plamb(np)*jumfac)
            lambda=max(lambda,plamb(np)/jumfac)
          else 
c
            open(48, file='chi.one', status='old')
            read(48,*) chione
c
            if(np.eq.1.and.pchi(1).lt.1.) then
              if(chione.eq.0) then
                lambda=lambda*redfac
              else
                power=1./float(2**chione)
                lambda=lambda*redfac**power
              end if
            else
              if(chione.eq.0) then
                lambda=lambda/redfac
              else
                power=1./float(2**chione)
                lambda=lambda/redfac**power
              end if
            end if
          end if
        end if
      end if
c
      rewind(14)
      write(14,75) lambda
75    format(e15.6)
c
      open(38, file='log.file',status='old',form='formatted')
155   read(38,55,end=98) alog
55    format(a1)
      go to 155
c
98    write(io,25) lambda
      backspace(38)
      write(38,25) lambda
25    format('value of lambda: ',e15.6)
c
      stop
      end
