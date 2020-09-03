c
c     version 1.0  Apr 1996
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |                   *****  I N I T 2  ****                     |
c     |                                                              |
c     |                  Initialize file lambda                      |
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
c        13 -- output: lambda
c
c
c     ----------------------------------------------------------------
c
c
c
      open(13, file='lambda')
      open(14, file='chi.one')
      open(15, file='stop.in')
c
      write(13,*) 0.
      write(14,*) 0
      write(15,*) 0
c
      stop
      end
