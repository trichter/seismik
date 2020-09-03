c
c     version 1.1  Apr 1995
c
c     2D stencils for FD
c
c     ----------------------------------------------------------------
c
      subroutine stencils2d(ixp,iyp,izp,ttrans)
c
      include 'fd.par'
      include 'fd.com'
c
      go to (1111,2222,9998,9998,5555,6666) iside
      go to 9998
c
1111  continue
c
c       stencils for left side of box
c
c          try stencil B
c
             if(max(time(ixp-1,iyp,izp+1),time(ixp,iyp,izp+1)).
     +       lt.1.e9) call timeb2d(ixp,iyp,izp+1,ixp-1,iyp,izp+1,
     +                             ixp,iyp,izp,ixp-1,iyp,izp,3)
             if(max(time(ixp-1,iyp,izp-1),time(ixp,iyp,izp-1)).
     +       lt.1.e9) call timeb2d(ixp,iyp,izp-1,ixp-1,iyp,izp-1,
     +                             ixp,iyp,izp,ixp-1,iyp,izp,3)
c
           if(tptmin.gt.1.e9) then      
c
c          try stencil C
c
           if(max(time(ixp,iyp,izp-1),time(ixp,iyp,izp+1)).
     +       lt.1.e9) 
     +          call timec2d(ixp,iyp,izp-1,ixp,iyp,izp,ixp,iyp,izp+1,
     +          ixp-1,iyp,izp,ixp-1,iyp,izp-1,ixp-1,iyp,izp+1)
           end if
c
c          try stencil C in 1D
c
           call time1d(ixp,iyp,izp,ixp-1,iyp,izp,6)
c
           ttrans=tptmin
c
c          try stencil C as a headwave
c
           if(time(ixp-1,iyp,izp+1).lt.1.e9) 
     +       call time1d(ixp-1,iyp,izp+1,ixp-1,iyp,izp,8)
           if(time(ixp-1,iyp,izp-1).lt.1.e9) 
     +       call time1d(ixp-1,iyp,izp-1,ixp-1,iyp,izp,8)
c
        return
c      
2222  continue
c
c       stencils for right side of box
c
c          try stencil B
c
             if(max(time(ixp+1,iyp,izp+1),time(ixp,iyp,izp+1)).
     +       lt.1.e9) call timeb2d(ixp,iyp,izp+1,ixp+1,iyp,izp+1,
     +                     ixp,iyp,izp,ixp+1,iyp,izp,3)
             if(max(time(ixp+1,iyp,izp-1),time(ixp,iyp,izp-1)).
     +       lt.1.e9) call timeb2d(ixp,iyp,izp-1,ixp+1,iyp,izp-1,
     +                             ixp,iyp,izp,ixp+1,iyp,izp,3)
c
           if(tptmin.gt.1.e9) then      
c
c          try stencil C
c
           if(max(time(ixp,iyp,izp-1),time(ixp,iyp,izp+1)).
     +       lt.1.e9) 
     +       call timec2d(ixp,iyp,izp-1,ixp,iyp,izp,ixp,iyp,izp+1,
     +          ixp+1,iyp,izp,ixp+1,iyp,izp-1,ixp+1,iyp,izp+1)
           end if
c
c          try stencil C in 1D
c
           call time1d(ixp,iyp,izp,ixp+1,iyp,izp,6)
c
           ttrans=tptmin
c
c          try stencil C as a headwave
c
           if(time(ixp+1,iyp,izp+1).lt.1.e9) 
     +       call time1d(ixp+1,iyp,izp+1,ixp+1,iyp,izp,8)
           if(time(ixp+1,iyp,izp-1).lt.1.e9) 
     +       call time1d(ixp+1,iyp,izp-1,ixp+1,iyp,izp,8)
c
        return
c      
5555  continue
c
c       stencils for top side of box
c
c          try stencil B
c
             if(max(time(ixp+1,iyp,izp-1),time(ixp+1,iyp,izp)).
     +       lt.1.e9) call timeb2d(ixp+1,iyp,izp,ixp+1,iyp,izp-1,
     +                             ixp,iyp,izp,ixp,iyp,izp-1,3)
             if(max(time(ixp-1,iyp,izp-1),time(ixp-1,iyp,izp)).
     +       lt.1.e9) call timeb2d(ixp-1,iyp,izp,ixp-1,iyp,izp-1,
     +                             ixp,iyp,izp,ixp,iyp,izp-1,3)
c
           if(tptmin.gt.1.e9) then      
c
c          try stencil C
c
           if(max(time(ixp-1,iyp,izp),time(ixp+1,iyp,izp)).
     +       lt.1.e9)           
     +       call timec2d(ixp-1,iyp,izp,ixp,iyp,izp,ixp+1,iyp,izp,
     +          ixp,iyp,izp-1,ixp-1,iyp,izp-1,ixp+1,iyp,izp-1)
           end if
c
c          try stencil C in 1D
c
           call time1d(ixp,iyp,izp,ixp,iyp,izp-1,6)
c
           ttrans=tptmin
c
c          try stencil C as a headwave
c
           if(time(ixp+1,iyp,izp-1).lt.1.e9) 
     +       call time1d(ixp+1,iyp,izp-1,ixp,iyp,izp-1,8)
           if(time(ixp-1,iyp,izp-1).lt.1.e9) 
     +       call time1d(ixp-1,iyp,izp-1,ixp,iyp,izp-1,8)
c
        return
c      
6666  continue
c
c       stencils for bottom side of box
c
c          try stencil B
c
             if(max(time(ixp+1,iyp,izp+1),time(ixp+1,iyp,izp)).
     +       lt.1.e9) call timeb2d(ixp+1,iyp,izp,ixp+1,iyp,izp+1,
     +                             ixp,iyp,izp,ixp,iyp,izp+1,3)
             if(max(time(ixp-1,iyp,izp+1),time(ixp-1,iyp,izp)).
     +       lt.1.e9) call timeb2d(ixp-1,iyp,izp,ixp-1,iyp,izp+1,
     +                             ixp,iyp,izp,ixp,iyp,izp+1,3)
c
           if(tptmin.gt.1.e9) then      
c
c          try stencil C
c
           if(max(time(ixp-1,iyp,izp),time(ixp+1,iyp,izp)).
     +       lt.1.e9)           
     +       call timec2d(ixp-1,iyp,izp,ixp,iyp,izp,ixp+1,iyp,izp,
     +          ixp,iyp,izp+1,ixp-1,iyp,izp+1,ixp+1,iyp,izp+1)
           end if
c
c          try stencil C in 1D
c
           call time1d(ixp,iyp,izp,ixp,iyp,izp+1,6)
c
           ttrans=tptmin
c
c          try stencil C as a headwave
c 
           if(time(ixp+1,iyp,izp+1).lt.1.e9) 
     +       call time1d(ixp+1,iyp,izp+1,ixp,iyp,izp+1,8)
           if(time(ixp-1,iyp,izp+1).lt.1.e9) 
     +       call time1d(ixp-1,iyp,izp+1,ixp,iyp,izp+1,8)
c
        return
c      
9998  write(6,95)
95    format(/'***  could not decide which side to expand  ***'/)
      stop   
      end
