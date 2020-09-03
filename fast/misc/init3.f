c
c     version 1.0  Apr 1996
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |                   *****  I N I T 3  ****                     |
c     |                                                              |
c     |  Initialize files current.iteration, chi.prev and time.out   |
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
c        11 -- output: current.iteration
c
c        12 -- output: chi.prev
c
c        14 -- output: time.out
c
c
c     ----------------------------------------------------------------
c
c
      open(11, file='current.iteration')
      open(12, file='chi.prev')
      open(13, file='stop.in')
      open(14, file='time.out')
c
      write(11,*) 1
      write(12,*) 1.e10
      write(13,*) 0
      write(14,*) 1
c
      stop
      end
