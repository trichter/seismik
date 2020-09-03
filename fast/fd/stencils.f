c 
c     version 1.1  Apr 1995
c
c     3D stencils for FD
c
c     ---------------------------------------------------------------- 
c
      subroutine stencils(ixp,iyp,izp,ttrans)
c
      include 'fd.par'
      include 'fd.com'
c
      go to (1111,2222,3333,4444,5555,6666) iside
c      
1111  continue
c
c       stencils for left side of box
c
c          try stencil A
c
           if(max(time(ixp-1,iyp+1,izp),time(ixp-1,iyp+1,izp+1),
     +        time(ixp-1,iyp,izp+1)).lt.tptmin) 
     +       call timea(ixp,iyp+1,izp+1,ixp-1,iyp,izp)
           if(max(time(ixp-1,iyp-1,izp),time(ixp-1,iyp-1,izp+1),
     +        time(ixp-1,iyp,izp+1)).lt.tptmin) 
     +       call timea(ixp,iyp-1,izp+1,ixp-1,iyp,izp)
           if(max(time(ixp-1,iyp+1,izp),time(ixp-1,iyp+1,izp-1),
     +        time(ixp-1,iyp,izp-1)).lt.tptmin) 
     +       call timea(ixp,iyp+1,izp-1,ixp-1,iyp,izp)
           if(max(time(ixp-1,iyp-1,izp),time(ixp-1,iyp-1,izp-1),
     +        time(ixp-1,iyp,izp-1)).lt.tptmin) 
     +       call timea(ixp,iyp-1,izp-1,ixp-1,iyp,izp)
c
           if(tptmin.gt.1.e9) then     
c
c          try stencil B
c
           if(max(time(ixp-1,iyp+1,izp),time(ixp,iyp+1,izp-1),
     +        time(ixp,iyp+1,izp+1)).lt.tptmin) 
     +       call timeb(ixp,iyp+1,izp-1,ixp,iyp+1,izp,ixp-1,iyp+1,izp,
     +          ixp,iyp+1,izp+1,ixp,iyp,izp,ixp-1,iyp,izp)
           if(max(time(ixp-1,iyp-1,izp),time(ixp,iyp-1,izp-1),
     +        time(ixp,iyp-1,izp+1)).lt.tptmin) 
     +       call timeb(ixp,iyp-1,izp-1,ixp,iyp-1,izp,ixp-1,iyp-1,izp,
     +          ixp,iyp-1,izp+1,ixp,iyp,izp,ixp-1,iyp,izp)
           if(max(time(ixp-1,iyp,izp+1),time(ixp,iyp-1,izp+1),
     +        time(ixp,iyp+1,izp+1)).lt.tptmin) 
     +       call timeb(ixp,iyp-1,izp+1,ixp,iyp,izp+1,ixp-1,iyp,izp+1,
     +          ixp,iyp+1,izp+1,ixp,iyp,izp,ixp-1,iyp,izp)
           if(max(time(ixp-1,iyp,izp-1),time(ixp,iyp-1,izp-1),
     +        time(ixp,iyp+1,izp-1)).lt.tptmin) 
     +       call timeb(ixp,iyp-1,izp-1,ixp,iyp,izp-1,ixp-1,iyp,izp-1,
     +          ixp,iyp+1,izp-1,ixp,iyp,izp,ixp-1,iyp,izp)
           end if
c
c          try stencil B in 2D
c 
           if(max(time(ixp-1,iyp,izp+1),time(ixp,iyp,izp+1)).lt.tptmin) 
     +     call timeb2d(ixp,iyp,izp+1,ixp-1,iyp,izp+1,
     +                  ixp,iyp,izp,ixp-1,iyp,izp,3)
           if(max(time(ixp-1,iyp,izp-1),time(ixp,iyp,izp-1)).lt.tptmin) 
     +     call timeb2d(ixp,iyp,izp-1,ixp-1,iyp,izp-1,
     +                  ixp,iyp,izp,ixp-1,iyp,izp,3)
           if(max(time(ixp-1,iyp+1,izp),time(ixp,iyp+1,izp)).lt.tptmin) 
     +     call timeb2d(ixp,iyp+1,izp,ixp-1,iyp+1,izp,
     +                  ixp,iyp,izp,ixp-1,iyp,izp,3)
           if(max(time(ixp-1,iyp-1,izp),time(ixp,iyp-1,izp)).lt.tptmin) 
     +     call timeb2d(ixp,iyp-1,izp,ixp-1,iyp-1,izp,
     +                  ixp,iyp,izp,ixp-1,iyp,izp,3)
c
           if(tptmin.gt.1.e9) then     
c
c          try stencil C
c
           ntype=0
           if(iyp.eq.side(3)) ntype=1
           if(iyp.eq.side(4)) ntype=2
           if(izp.eq.side(5)) ntype=3
           if(izp.eq.side(6)) ntype=4
c
           if(ntype.eq.0) 
     +       call timec(ixp,iyp-1,izp,ixp,iyp,izp-1,ixp,iyp,izp,
     +          ixp,iyp,izp+1,ixp,iyp+1,izp,ixp-1,iyp,izp,
     +          ixp,iyp-1,izp-1,ixp,iyp-1,izp+1,
     +          ixp,iyp+1,izp-1,ixp,iyp+1,izp+1)
           end if
c
c          try stencil C in 1D
c
           call time1d(ixp,iyp,izp,ixp-1,iyp,izp,6)
c
           ttrans=tptmin
c 
c          try stencil B in 2D as a headwave
c
           if(max(time(ixp-1,iyp,izp+1),time(ixp-1,iyp+1,izp),
     +            time(ixp-1,iyp+1,izp+1)).lt.tptmin) 
     +       call timeb2d(ixp-1,iyp+1,izp+1,ixp-1,iyp+1,izp,
     +                  ixp-1,iyp,izp+1,ixp-1,iyp,izp,7)
           if(max(time(ixp-1,iyp,izp-1),time(ixp-1,iyp+1,izp),
     +            time(ixp-1,iyp+1,izp-1)).lt.tptmin) 
     +       call timeb2d(ixp-1,iyp+1,izp-1,ixp-1,iyp+1,izp,
     +                  ixp-1,iyp,izp-1,ixp-1,iyp,izp,7)
           if(max(time(ixp-1,iyp,izp+1),time(ixp-1,iyp-1,izp),
     +            time(ixp-1,iyp-1,izp+1)).lt.tptmin) 
     +       call timeb2d(ixp-1,iyp-1,izp+1,ixp-1,iyp-1,izp,
     +                  ixp-1,iyp,izp+1,ixp-1,iyp,izp,7)
           if(max(time(ixp-1,iyp,izp-1),time(ixp-1,iyp-1,izp),
     +            time(ixp-1,iyp-1,izp-1)).lt.tptmin) 
     +       call timeb2d(ixp-1,iyp-1,izp-1,ixp-1,iyp-1,izp,
     +                  ixp-1,iyp,izp-1,ixp-1,iyp,izp,7)
c
c          try stencil C in 1D as a headwave
c
           if(time(ixp-1,iyp,izp+1).lt.tptmin) 
     +       call time1d(ixp-1,iyp,izp+1,ixp-1,iyp,izp,8)
           if(time(ixp-1,iyp,izp-1).lt.tptmin) 
     +       call time1d(ixp-1,iyp,izp-1,ixp-1,iyp,izp,8)
           if(time(ixp-1,iyp+1,izp).lt.tptmin) 
     +       call time1d(ixp-1,iyp+1,izp,ixp-1,iyp,izp,8)
           if(time(ixp-1,iyp-1,izp).lt.tptmin) 
     +       call time1d(ixp-1,iyp-1,izp,ixp-1,iyp,izp,8)
c
        return      
c      
2222  continue
c
c       stencils for right side of box
c
c          try stencil A
c
             if(max(time(ixp+1,iyp+1,izp),time(ixp+1,iyp+1,izp+1),
     +          time(ixp+1,iyp,izp+1)).lt.tptmin) 
     +         call timea(ixp,iyp+1,izp+1,ixp+1,iyp,izp)
             if(max(time(ixp+1,iyp-1,izp),time(ixp+1,iyp-1,izp+1),
     +          time(ixp+1,iyp,izp+1)).lt.tptmin) 
     +         call timea(ixp,iyp-1,izp+1,ixp+1,iyp,izp)
             if(max(time(ixp+1,iyp+1,izp),time(ixp+1,iyp+1,izp-1),
     +          time(ixp+1,iyp,izp-1)).lt.tptmin) 
     +         call timea(ixp,iyp+1,izp-1,ixp+1,iyp,izp)
             if(max(time(ixp+1,iyp-1,izp),time(ixp+1,iyp-1,izp-1),
     +          time(ixp+1,iyp,izp-1)).lt.tptmin) 
     +         call timea(ixp,iyp-1,izp-1,ixp+1,iyp,izp)
c
           if(tptmin.gt.1.e9) then     
c
c          try stencil B
c
             if(max(time(ixp+1,iyp+1,izp),time(ixp,iyp+1,izp-1),
     +          time(ixp,iyp+1,izp+1)).lt.tptmin) 
     +         call timeb(ixp,iyp+1,izp-1,ixp,iyp+1,izp,ixp+1,iyp+1,izp,
     +            ixp,iyp+1,izp+1,ixp,iyp,izp,ixp+1,iyp,izp)
             if(max(time(ixp+1,iyp-1,izp),time(ixp,iyp-1,izp-1),
     +          time(ixp,iyp-1,izp+1)).lt.tptmin) 
     +         call timeb(ixp,iyp-1,izp-1,ixp,iyp-1,izp,ixp+1,iyp-1,izp,
     +            ixp,iyp-1,izp+1,ixp,iyp,izp,ixp+1,iyp,izp)
             if(max(time(ixp+1,iyp,izp+1),time(ixp,iyp-1,izp+1),
     +          time(ixp,iyp+1,izp+1)).lt.tptmin) 
     +         call timeb(ixp,iyp-1,izp+1,ixp,iyp,izp+1,ixp+1,iyp,izp+1,
     +            ixp,iyp+1,izp+1,ixp,iyp,izp,ixp+1,iyp,izp)
             if(max(time(ixp+1,iyp,izp-1),time(ixp,iyp-1,izp-1),
     +          time(ixp,iyp+1,izp-1)).lt.tptmin) 
     +         call timeb(ixp,iyp-1,izp-1,ixp,iyp,izp-1,ixp+1,iyp,izp-1,
     +            ixp,iyp+1,izp-1,ixp,iyp,izp,ixp+1,iyp,izp)
           end if
c
c          try stencil B in 2D
c
           if(max(time(ixp+1,iyp,izp+1),time(ixp,iyp,izp+1)).lt.tptmin) 
     +     call timeb2d(ixp,iyp,izp+1,ixp+1,iyp,izp+1,
     +                  ixp,iyp,izp,ixp+1,iyp,izp,3)
           if(max(time(ixp+1,iyp,izp-1),time(ixp,iyp,izp-1)).lt.tptmin) 
     +     call timeb2d(ixp,iyp,izp-1,ixp+1,iyp,izp-1,
     +                  ixp,iyp,izp,ixp+1,iyp,izp,3)
           if(max(time(ixp+1,iyp+1,izp),time(ixp,iyp+1,izp)).lt.tptmin) 
     +     call timeb2d(ixp,iyp+1,izp,ixp+1,iyp+1,izp,
     +                  ixp,iyp,izp,ixp+1,iyp,izp,3)
           if(max(time(ixp+1,iyp-1,izp),time(ixp,iyp-1,izp)).lt.tptmin) 
     +     call timeb2d(ixp,iyp-1,izp,ixp+1,iyp-1,izp,
     +                  ixp,iyp,izp,ixp+1,iyp,izp,3)
c
           if(tptmin.gt.1.e9) then     
c
c          try stencil C
c
           ntype=0
           if(iyp.eq.side(3)) ntype=1
           if(iyp.eq.side(4)) ntype=2
           if(izp.eq.side(5)) ntype=3
           if(izp.eq.side(6)) ntype=4
c
           if(ntype.eq.0) 
     +       call timec(ixp,iyp-1,izp,ixp,iyp,izp-1,ixp,iyp,izp,
     +          ixp,iyp,izp+1,ixp,iyp+1,izp,ixp+1,iyp,izp,
     +          ixp,iyp-1,izp-1,ixp,iyp-1,izp+1,
     +          ixp,iyp+1,izp-1,ixp,iyp+1,izp+1)
           end if
c
c          try stencil C in 1D
c
           call time1d(ixp,iyp,izp,ixp+1,iyp,izp,6)
c
           ttrans=tptmin
c
c          try stencil B in 2D as a headwave
c
           if(max(time(ixp+1,iyp,izp+1),time(ixp+1,iyp+1,izp),
     +            time(ixp+1,iyp+1,izp+1)).lt.tptmin) 
     +       call timeb2d(ixp+1,iyp+1,izp+1,ixp+1,iyp+1,izp,
     +                  ixp+1,iyp,izp+1,ixp+1,iyp,izp,7)
           if(max(time(ixp+1,iyp,izp-1),time(ixp+1,iyp+1,izp),
     +            time(ixp+1,iyp+1,izp-1)).lt.tptmin) 
     +       call timeb2d(ixp+1,iyp+1,izp-1,ixp+1,iyp+1,izp,
     +                  ixp+1,iyp,izp-1,ixp+1,iyp,izp,7)
           if(max(time(ixp+1,iyp,izp+1),time(ixp+1,iyp-1,izp),
     +            time(ixp+1,iyp-1,izp+1)).lt.tptmin) 
     +       call timeb2d(ixp+1,iyp-1,izp+1,ixp+1,iyp-1,izp,
     +                  ixp+1,iyp,izp+1,ixp+1,iyp,izp,7)
           if(max(time(ixp+1,iyp,izp-1),time(ixp+1,iyp-1,izp),
     +            time(ixp+1,iyp-1,izp-1)).lt.tptmin) 
     +       call timeb2d(ixp+1,iyp-1,izp-1,ixp+1,iyp-1,izp,
     +                  ixp+1,iyp,izp-1,ixp+1,iyp,izp,7)
c
c          try stencil C in 1D as a headwave
c
           if(time(ixp+1,iyp,izp+1).lt.tptmin) 
     +       call time1d(ixp+1,iyp,izp+1,ixp+1,iyp,izp,8)
           if(time(ixp+1,iyp,izp-1).lt.tptmin) 
     +       call time1d(ixp+1,iyp,izp-1,ixp+1,iyp,izp,8)
           if(time(ixp+1,iyp+1,izp).lt.tptmin) 
     +       call time1d(ixp+1,iyp+1,izp,ixp+1,iyp,izp,8)
           if(time(ixp+1,iyp-1,izp).lt.tptmin) 
     +       call time1d(ixp+1,iyp-1,izp,ixp+1,iyp,izp,8)
c
        return      
c      
3333  continue
c
c       stencils for back side of box
c
c          try stencil A
c
             if(max(time(ixp+1,iyp-1,izp),time(ixp+1,iyp-1,izp+1),
     +          time(ixp,iyp-1,izp+1)).lt.tptmin) 
     +         call timea(ixp+1,iyp,izp+1,ixp,iyp-1,izp)
             if(max(time(ixp-1,iyp-1,izp),time(ixp-1,iyp-1,izp+1),
     +          time(ixp,iyp-1,izp+1)).lt.tptmin) 
     +         call timea(ixp-1,iyp,izp+1,ixp,iyp-1,izp)
             if(max(time(ixp+1,iyp-1,izp),time(ixp+1,iyp-1,izp-1),
     +          time(ixp,iyp-1,izp-1)).lt.tptmin) 
     +         call timea(ixp+1,iyp,izp-1,ixp,iyp-1,izp)
             if(max(time(ixp-1,iyp-1,izp),time(ixp-1,iyp-1,izp-1),
     +          time(ixp,iyp-1,izp-1)).lt.tptmin) 
     +         call timea(ixp-1,iyp,izp-1,ixp,iyp-1,izp)
c
           if(tptmin.gt.1.e9) then     
c
c          try stencil B
c
             if(max(time(ixp+1,iyp-1,izp),time(ixp+1,iyp,izp-1),
     +          time(ixp+1,iyp,izp+1)).lt.tptmin) 
     +         call timeb(ixp+1,iyp,izp-1,ixp+1,iyp,izp,ixp+1,iyp-1,izp,
     +            ixp+1,iyp,izp+1,ixp,iyp,izp,ixp,iyp-1,izp)
             if(max(time(ixp-1,iyp-1,izp),time(ixp-1,iyp,izp-1),
     +          time(ixp-1,iyp,izp+1)).lt.tptmin) 
     +         call timeb(ixp-1,iyp,izp-1,ixp-1,iyp,izp,ixp-1,iyp-1,izp,
     +            ixp-1,iyp,izp+1,ixp,iyp,izp,ixp,iyp-1,izp)
             if(max(time(ixp,iyp-1,izp+1),time(ixp-1,iyp,izp+1),
     +          time(ixp+1,iyp,izp+1)).lt.tptmin) 
     +         call timeb(ixp-1,iyp,izp+1,ixp,iyp,izp+1,ixp,iyp-1,izp+1,
     +            ixp+1,iyp,izp+1,ixp,iyp,izp,ixp,iyp-1,izp)
             if(max(time(ixp,iyp-1,izp-1),time(ixp-1,iyp,izp-1),
     +          time(ixp+1,iyp,izp-1)).lt.tptmin) 
     +         call timeb(ixp-1,iyp,izp-1,ixp,iyp,izp-1,ixp,iyp-1,izp-1,
     +            ixp+1,iyp,izp-1,ixp,iyp,izp,ixp,iyp-1,izp)
           end if
c
c          try stencil B in 2D
c
           if(max(time(ixp,iyp-1,izp+1),time(ixp,iyp,izp+1)).lt.tptmin) 
     +     call timeb2d(ixp,iyp,izp+1,ixp,iyp-1,izp+1,
     +                  ixp,iyp,izp,ixp,iyp-1,izp,3)
           if(max(time(ixp,iyp-1,izp-1),time(ixp,iyp,izp-1)).lt.tptmin) 
     +     call timeb2d(ixp,iyp,izp-1,ixp,iyp-1,izp-1,
     +                  ixp,iyp,izp,ixp,iyp-1,izp,3)
           if(max(time(ixp+1,iyp-1,izp),time(ixp+1,iyp,izp)).lt.tptmin) 
     +     call timeb2d(ixp+1,iyp,izp,ixp+1,iyp-1,izp,
     +                  ixp,iyp,izp,ixp,iyp-1,izp,3)
           if(max(time(ixp-1,iyp-1,izp),time(ixp-1,iyp,izp)).lt.tptmin) 
     +     call timeb2d(ixp-1,iyp,izp,ixp-1,iyp-1,izp,
     +                  ixp,iyp,izp,ixp,iyp-1,izp,3)
c
           if(tptmin.gt.1.e9) then     
c
c          try stencil C
c
           ntype=0
           if(ixp.eq.side(1)) ntype=1
           if(ixp.eq.side(2)) ntype=2
           if(izp.eq.side(5)) ntype=3
           if(izp.eq.side(6)) ntype=4
c
           if(ntype.eq.0)
     +       call timec(ixp-1,iyp,izp,ixp,iyp,izp-1,ixp,iyp,izp,
     +          ixp,iyp,izp+1,ixp+1,iyp,izp,ixp,iyp-1,izp,
     +          ixp-1,iyp,izp-1,ixp-1,iyp,izp+1,
     +          ixp+1,iyp,izp-1,ixp+1,iyp,izp+1)
           end if
c
c          try stencil C in 1D
c
           call time1d(ixp,iyp,izp,ixp,iyp-1,izp,6)
c
           ttrans=tptmin
c
c          try stencil B in 2D as a headwave
c
           if(max(time(ixp,iyp-1,izp+1),time(ixp+1,iyp-1,izp),
     +            time(ixp+1,iyp-1,izp+1)).lt.tptmin) 
     +       call timeb2d(ixp+1,iyp-1,izp+1,ixp+1,iyp-1,izp,
     +                  ixp,iyp-1,izp+1,ixp,iyp-1,izp,7)
           if(max(time(ixp,iyp-1,izp-1),time(ixp+1,iyp-1,izp),
     +            time(ixp+1,iyp-1,izp-1)).lt.tptmin) 
     +       call timeb2d(ixp+1,iyp-1,izp-1,ixp+1,iyp-1,izp,
     +                  ixp,iyp-1,izp-1,ixp,iyp-1,izp,7)
           if(max(time(ixp,iyp-1,izp+1),time(ixp-1,iyp-1,izp),
     +            time(ixp-1,iyp-1,izp+1)).lt.tptmin) 
     +       call timeb2d(ixp-1,iyp-1,izp+1,ixp-1,iyp-1,izp,
     +                  ixp,iyp-1,izp+1,ixp,iyp-1,izp,7)
           if(max(time(ixp,iyp-1,izp-1),time(ixp-1,iyp-1,izp),
     +            time(ixp-1,iyp-1,izp-1)).lt.tptmin) 
     +       call timeb2d(ixp-1,iyp-1,izp-1,ixp-1,iyp-1,izp,
     +                  ixp,iyp-1,izp-1,ixp,iyp-1,izp,7)
c
c          try stencil C in 1D as a headwave
c
           if(time(ixp,iyp-1,izp+1).lt.tptmin) 
     +       call time1d(ixp,iyp-1,izp+1,ixp,iyp-1,izp,8)
           if(time(ixp,iyp-1,izp-1).lt.tptmin) 
     +       call time1d(ixp,iyp-1,izp-1,ixp,iyp-1,izp,8)
           if(time(ixp+1,iyp-1,izp).lt.tptmin) 
     +       call time1d(ixp+1,iyp-1,izp,ixp,iyp-1,izp,8)
           if(time(ixp-1,iyp-1,izp).lt.tptmin) 
     +       call time1d(ixp-1,iyp-1,izp,ixp,iyp-1,izp,8)
c
        return      
c
4444  continue
c
c       stencils for front side of box
c
c          try stencil A
c
             if(max(time(ixp+1,iyp+1,izp),time(ixp+1,iyp+1,izp+1),
     +          time(ixp,iyp+1,izp+1)).lt.tptmin) 
     +         call timea(ixp+1,iyp,izp+1,ixp,iyp+1,izp)
             if(max(time(ixp-1,iyp+1,izp),time(ixp-1,iyp+1,izp+1),
     +          time(ixp,iyp+1,izp+1)).lt.tptmin) 
     +         call timea(ixp-1,iyp,izp+1,ixp,iyp+1,izp)
             if(max(time(ixp+1,iyp+1,izp),time(ixp+1,iyp+1,izp-1),
     +          time(ixp,iyp+1,izp-1)).lt.tptmin) 
     +         call timea(ixp+1,iyp,izp-1,ixp,iyp+1,izp)
             if(max(time(ixp-1,iyp+1,izp),time(ixp-1,iyp+1,izp-1),
     +          time(ixp,iyp+1,izp-1)).lt.tptmin) 
     +         call timea(ixp-1,iyp,izp-1,ixp,iyp+1,izp)
c
           if(tptmin.gt.1.e9) then     
c
c          try stencil B
c
             if(max(time(ixp+1,iyp+1,izp),time(ixp+1,iyp,izp-1),
     +          time(ixp+1,iyp,izp+1)).lt.tptmin) 
     +         call timeb(ixp+1,iyp,izp-1,ixp+1,iyp,izp,ixp+1,iyp+1,izp,
     +            ixp+1,iyp,izp+1,ixp,iyp,izp,ixp,iyp+1,izp)
             if(max(time(ixp-1,iyp+1,izp),time(ixp-1,iyp,izp-1),
     +          time(ixp-1,iyp,izp+1)).lt.tptmin) 
     +         call timeb(ixp-1,iyp,izp-1,ixp-1,iyp,izp,ixp-1,iyp+1,izp,
     +            ixp-1,iyp,izp+1,ixp,iyp,izp,ixp,iyp+1,izp)
             if(max(time(ixp,iyp+1,izp+1),time(ixp-1,iyp,izp+1),
     +          time(ixp+1,iyp,izp+1)).lt.tptmin) 
     +         call timeb(ixp-1,iyp,izp+1,ixp,iyp,izp+1,ixp,iyp+1,izp+1,
     +            ixp+1,iyp,izp+1,ixp,iyp,izp,ixp,iyp+1,izp)
             if(max(time(ixp,iyp+1,izp-1),time(ixp-1,iyp,izp-1),
     +          time(ixp+1,iyp,izp-1)).lt.tptmin) 
     +         call timeb(ixp-1,iyp,izp-1,ixp,iyp,izp-1,ixp,iyp+1,izp-1,
     +            ixp+1,iyp,izp-1,ixp,iyp,izp,ixp,iyp+1,izp)
           end if
c
c          try stencil B in 2D
c
           if(max(time(ixp,iyp+1,izp+1),time(ixp,iyp,izp+1)).lt.tptmin) 
     +     call timeb2d(ixp,iyp,izp+1,ixp,iyp+1,izp+1,
     +                  ixp,iyp,izp,ixp,iyp+1,izp,3)
           if(max(time(ixp,iyp+1,izp-1),time(ixp,iyp,izp-1)).lt.tptmin) 
     +     call timeb2d(ixp,iyp,izp-1,ixp,iyp+1,izp-1,
     +                  ixp,iyp,izp,ixp,iyp+1,izp,3)
           if(max(time(ixp+1,iyp+1,izp),time(ixp+1,iyp,izp)).lt.tptmin) 
     +     call timeb2d(ixp+1,iyp,izp,ixp+1,iyp+1,izp,
     +                  ixp,iyp,izp,ixp,iyp+1,izp,3)
           if(max(time(ixp-1,iyp+1,izp),time(ixp-1,iyp,izp)).lt.tptmin) 
     +     call timeb2d(ixp-1,iyp,izp,ixp-1,iyp+1,izp,
     +                  ixp,iyp,izp,ixp,iyp+1,izp,3)
c
           if(tptmin.gt.1.e9) then     
c
c          try stencil C
c
           ntype=0
           if(ixp.eq.side(1)) ntype=1
           if(ixp.eq.side(2)) ntype=2
           if(izp.eq.side(5)) ntype=3
           if(izp.eq.side(6)) ntype=4
c
           if(ntype.eq.0) 
     +       call timec(ixp-1,iyp,izp,ixp,iyp,izp-1,ixp,iyp,izp,
     +          ixp,iyp,izp+1,ixp+1,iyp,izp,ixp,iyp+1,izp,
     +          ixp-1,iyp,izp-1,ixp-1,iyp,izp+1,
     +          ixp+1,iyp,izp-1,ixp+1,iyp,izp+1)
           end if
c
c          try stencil C in 1D
c
           call time1d(ixp,iyp,izp,ixp,iyp+1,izp,6)
c
           ttrans=tptmin
c
c          try stencil B in 2D as a headwave
c
           if(max(time(ixp,iyp+1,izp+1),time(ixp+1,iyp+1,izp),
     +            time(ixp+1,iyp+1,izp+1)).lt.tptmin) 
     +       call timeb2d(ixp+1,iyp+1,izp+1,ixp+1,iyp+1,izp,
     +                  ixp,iyp+1,izp+1,ixp,iyp+1,izp,7)
           if(max(time(ixp,iyp+1,izp-1),time(ixp+1,iyp+1,izp),
     +            time(ixp+1,iyp+1,izp-1)).lt.tptmin) 
     +       call timeb2d(ixp+1,iyp+1,izp-1,ixp+1,iyp+1,izp,
     +                  ixp,iyp+1,izp-1,ixp,iyp+1,izp,7)
           if(max(time(ixp,iyp+1,izp+1),time(ixp-1,iyp+1,izp),
     +            time(ixp-1,iyp+1,izp+1)).lt.tptmin) 
     +       call timeb2d(ixp-1,iyp+1,izp+1,ixp-1,iyp+1,izp,
     +                  ixp,iyp+1,izp+1,ixp,iyp+1,izp,7)
           if(max(time(ixp,iyp+1,izp-1),time(ixp-1,iyp+1,izp),
     +            time(ixp-1,iyp+1,izp-1)).lt.tptmin) 
     +       call timeb2d(ixp-1,iyp+1,izp-1,ixp-1,iyp+1,izp,
     +                  ixp,iyp+1,izp-1,ixp,iyp+1,izp,7)
c
c          try stencil C in 1D as a headwave
c
           if(time(ixp,iyp+1,izp+1).lt.tptmin) 
     +       call time1d(ixp,iyp+1,izp+1,ixp,iyp+1,izp,8)
           if(time(ixp,iyp+1,izp-1).lt.tptmin) 
     +       call time1d(ixp,iyp+1,izp-1,ixp,iyp+1,izp,8)
           if(time(ixp+1,iyp+1,izp).lt.tptmin) 
     +       call time1d(ixp+1,iyp+1,izp,ixp,iyp+1,izp,8)
           if(time(ixp-1,iyp+1,izp).lt.tptmin) 
     +       call time1d(ixp-1,iyp+1,izp,ixp,iyp+1,izp,8)
c
        return      
c      
5555  continue
c
c       stencils for top side of box
c
c          try stencil A
c
             if(max(time(ixp,iyp+1,izp-1),time(ixp+1,iyp+1,izp-1),
     +          time(ixp+1,iyp,izp-1)).lt.tptmin) 
     +         call timea(ixp+1,iyp+1,izp,ixp,iyp,izp-1)
             if(max(time(ixp,iyp-1,izp-1),time(ixp+1,iyp-1,izp-1),
     +          time(ixp+1,iyp,izp-1)).lt.tptmin) 
     +         call timea(ixp+1,iyp-1,izp,ixp,iyp,izp-1)
             if(max(time(ixp,iyp+1,izp-1),time(ixp-1,iyp+1,izp-1),
     +          time(ixp-1,iyp,izp-1)).lt.tptmin) 
     +         call timea(ixp-1,iyp+1,izp,ixp,iyp,izp-1)
             if(max(time(ixp,iyp-1,izp-1),time(ixp-1,iyp-1,izp-1),
     +          time(ixp-1,iyp,izp-1)).lt.tptmin) 
     +         call timea(ixp-1,iyp-1,izp,ixp,iyp,izp-1)
c
           if(tptmin.gt.1.e9) then     
c
c          try stencil B
c
             if(max(time(ixp,iyp+1,izp-1),time(ixp-1,iyp+1,izp),
     +          time(ixp+1,iyp+1,izp)).lt.tptmin) 
     +         call timeb(ixp-1,iyp+1,izp,ixp,iyp+1,izp,ixp,iyp+1,izp-1,
     +            ixp+1,iyp+1,izp,ixp,iyp,izp,ixp,iyp,izp-1)
             if(max(time(ixp,iyp-1,izp-1),time(ixp-1,iyp-1,izp),
     +          time(ixp+1,iyp-1,izp)).lt.tptmin) 
     +         call timeb(ixp-1,iyp-1,izp,ixp,iyp-1,izp,ixp,iyp-1,izp-1,
     +            ixp+1,iyp-1,izp,ixp,iyp,izp,ixp,iyp,izp-1)
             if(max(time(ixp+1,iyp,izp-1),time(ixp+1,iyp-1,izp),
     +          time(ixp+1,iyp+1,izp)).lt.tptmin) 
     +         call timeb(ixp+1,iyp-1,izp,ixp+1,iyp,izp,ixp+1,iyp,izp-1,
     +            ixp+1,iyp+1,izp,ixp,iyp,izp,ixp,iyp,izp-1)
             if(max(time(ixp-1,iyp,izp-1),time(ixp-1,iyp-1,izp),
     +          time(ixp-1,iyp+1,izp)).lt.tptmin) 
     +         call timeb(ixp-1,iyp-1,izp,ixp-1,iyp,izp,ixp-1,iyp,izp-1,
     +            ixp-1,iyp+1,izp,ixp,iyp,izp,ixp,iyp,izp-1)
           end if
c
c          try stencil B in 2D
c
           if(max(time(ixp+1,iyp,izp-1),time(ixp+1,iyp,izp)).lt.tptmin) 
     +     call timeb2d(ixp+1,iyp,izp,ixp+1,iyp,izp-1,
     +                  ixp,iyp,izp,ixp,iyp,izp-1,3)
           if(max(time(ixp-1,iyp,izp-1),time(ixp-1,iyp,izp)).lt.tptmin) 
     +     call timeb2d(ixp-1,iyp,izp,ixp-1,iyp,izp-1,
     +                  ixp,iyp,izp,ixp,iyp,izp-1,3)
           if(max(time(ixp,iyp+1,izp-1),time(ixp,iyp+1,izp)).lt.tptmin) 
     +     call timeb2d(ixp,iyp+1,izp,ixp,iyp+1,izp-1,
     +                  ixp,iyp,izp,ixp,iyp,izp-1,3)
           if(max(time(ixp,iyp-1,izp-1),time(ixp,iyp-1,izp)).lt.tptmin) 
     +     call timeb2d(ixp,iyp-1,izp,ixp,iyp-1,izp-1,
     +                  ixp,iyp,izp,ixp,iyp,izp-1,3)
c
           if(tptmin.gt.1.e9) then     
c
c          try stencil C
c
           ntype=0
           if(ixp.eq.side(1)) ntype=1
           if(ixp.eq.side(2)) ntype=2
           if(iyp.eq.side(3)) ntype=3
           if(iyp.eq.side(4)) ntype=4
c
           if(ntype.eq.0) 
     +       call timec(ixp,iyp-1,izp,ixp-1,iyp,izp,ixp,iyp,izp,
     +          ixp+1,iyp,izp,ixp,iyp+1,izp,ixp,iyp,izp-1,
     +          ixp-1,iyp-1,izp,ixp-1,iyp+1,izp,
     +          ixp+1,iyp-1,izp,ixp+1,iyp+1,izp)
           end if
c
c          try stencil C in 1D
c
           call time1d(ixp,iyp,izp,ixp,iyp,izp-1,6)
c
           ttrans=tptmin
c
c          try stencil B in 2D as a headwave
c
           if(max(time(ixp+1,iyp,izp-1),time(ixp,iyp+1,izp-1),
     +            time(ixp+1,iyp+1,izp-1)).lt.tptmin) 
     +       call timeb2d(ixp+1,iyp+1,izp-1,ixp,iyp+1,izp-1,
     +                  ixp+1,iyp,izp-1,ixp,iyp,izp-1,7)
           if(max(time(ixp-1,iyp,izp-1),time(ixp,iyp+1,izp-1),
     +            time(ixp-1,iyp+1,izp-1)).lt.tptmin) 
     +       call timeb2d(ixp-1,iyp+1,izp-1,ixp,iyp+1,izp-1,
     +                  ixp-1,iyp,izp-1,ixp,iyp,izp-1,7)
           if(max(time(ixp+1,iyp,izp-1),time(ixp,iyp-1,izp-1),
     +            time(ixp+1,iyp-1,izp-1)).lt.tptmin) 
     +       call timeb2d(ixp+1,iyp-1,izp-1,ixp,iyp-1,izp-1,
     +                  ixp+1,iyp,izp-1,ixp,iyp,izp-1,7)
           if(max(time(ixp-1,iyp,izp-1),time(ixp,iyp-1,izp-1),
     +            time(ixp-1,iyp-1,izp-1)).lt.tptmin) 
     +       call timeb2d(ixp-1,iyp-1,izp-1,ixp,iyp-1,izp-1,
     +                  ixp-1,iyp,izp-1,ixp,iyp,izp-1,7)
c
c          try stencil C in 1D as a headwave
c
           if(time(ixp+1,iyp,izp-1).lt.tptmin) 
     +       call time1d(ixp+1,iyp,izp-1,ixp,iyp,izp-1,8)
           if(time(ixp-1,iyp,izp-1).lt.tptmin) 
     +       call time1d(ixp-1,iyp,izp-1,ixp,iyp,izp-1,8)
           if(time(ixp,iyp+1,izp-1).lt.tptmin) 
     +       call time1d(ixp,iyp+1,izp-1,ixp,iyp,izp-1,8)
           if(time(ixp,iyp-1,izp-1).lt.tptmin) 
     +       call time1d(ixp,iyp-1,izp-1,ixp,iyp,izp-1,8)
c
        return      
c      
6666  continue
c
c       stencils for bottom side of box
c
c          try stencil A
c
             if(max(time(ixp,iyp+1,izp+1),time(ixp+1,iyp+1,izp+1),
     +          time(ixp+1,iyp,izp+1)).lt.tptmin) 
     +         call timea(ixp+1,iyp+1,izp,ixp,iyp,izp+1)
             if(max(time(ixp,iyp-1,izp+1),time(ixp+1,iyp-1,izp+1),
     +          time(ixp+1,iyp,izp+1)).lt.tptmin) 
     +         call timea(ixp+1,iyp-1,izp,ixp,iyp,izp+1)
             if(max(time(ixp,iyp+1,izp+1),time(ixp-1,iyp+1,izp+1),
     +          time(ixp-1,iyp,izp+1)).lt.tptmin) 
     +         call timea(ixp-1,iyp+1,izp,ixp,iyp,izp+1)
             if(max(time(ixp,iyp-1,izp+1),time(ixp-1,iyp-1,izp+1),
     +          time(ixp-1,iyp,izp+1)).lt.tptmin) 
     +         call timea(ixp-1,iyp-1,izp,ixp,iyp,izp+1)
c
           if(tptmin.gt.1.e9) then     
c
c          try stencil B
c
             if(max(time(ixp,iyp+1,izp+1),time(ixp-1,iyp+1,izp),
     +          time(ixp+1,iyp+1,izp)).lt.tptmin) 
     +         call timeb(ixp-1,iyp+1,izp,ixp,iyp+1,izp,ixp,iyp+1,izp+1,
     +            ixp+1,iyp+1,izp,ixp,iyp,izp,ixp,iyp,izp+1)
             if(max(time(ixp,iyp-1,izp+1),time(ixp-1,iyp-1,izp),
     +          time(ixp+1,iyp-1,izp)).lt.tptmin) 
     +         call timeb(ixp-1,iyp-1,izp,ixp,iyp-1,izp,ixp,iyp-1,izp+1,
     +            ixp+1,iyp-1,izp,ixp,iyp,izp,ixp,iyp,izp+1)
             if(max(time(ixp+1,iyp,izp+1),time(ixp+1,iyp-1,izp),
     +          time(ixp+1,iyp+1,izp)).lt.tptmin) 
     +         call timeb(ixp+1,iyp-1,izp,ixp+1,iyp,izp,ixp+1,iyp,izp+1,
     +            ixp+1,iyp+1,izp,ixp,iyp,izp,ixp,iyp,izp+1)
             if(max(time(ixp-1,iyp,izp+1),time(ixp-1,iyp-1,izp),
     +          time(ixp-1,iyp+1,izp)).lt.tptmin) 
     +         call timeb(ixp-1,iyp-1,izp,ixp-1,iyp,izp,ixp-1,iyp,izp+1,
     +            ixp-1,iyp+1,izp,ixp,iyp,izp,ixp,iyp,izp+1)
           end if
c
c          try stencil B in 2D
c
           if(max(time(ixp+1,iyp,izp+1),time(ixp+1,iyp,izp)).lt.tptmin) 
     +     call timeb2d(ixp+1,iyp,izp,ixp+1,iyp,izp+1,
     +                  ixp,iyp,izp,ixp,iyp,izp+1,3)
           if(max(time(ixp-1,iyp,izp+1),time(ixp-1,iyp,izp)).lt.tptmin) 
     +     call timeb2d(ixp-1,iyp,izp,ixp-1,iyp,izp+1,
     +                  ixp,iyp,izp,ixp,iyp,izp+1,3)
           if(max(time(ixp,iyp+1,izp+1),time(ixp,iyp+1,izp)).lt.tptmin) 
     +     call timeb2d(ixp,iyp+1,izp,ixp,iyp+1,izp+1,
     +                  ixp,iyp,izp,ixp,iyp,izp+1,3)
           if(max(time(ixp,iyp-1,izp+1),time(ixp,iyp-1,izp)).lt.tptmin) 
     +     call timeb2d(ixp,iyp-1,izp,ixp,iyp-1,izp+1,
     +                  ixp,iyp,izp,ixp,iyp,izp+1,3)
c
           if(tptmin.gt.1.e9) then     
c
c          try stencil C
c
           ntype=0
           if(ixp.eq.side(1)) ntype=1
           if(ixp.eq.side(2)) ntype=2
           if(iyp.eq.side(3)) ntype=3
           if(iyp.eq.side(4)) ntype=4
c
           if(ntype.eq.0) 
     +       call timec(ixp,iyp-1,izp,ixp-1,iyp,izp,ixp,iyp,izp,
     +          ixp+1,iyp,izp,ixp,iyp+1,izp,ixp,iyp,izp+1,
     +          ixp-1,iyp-1,izp,ixp-1,iyp+1,izp,
     +          ixp+1,iyp-1,izp,ixp+1,iyp+1,izp)
           end if
c
c          try stencil C in 1D
c
           call time1d(ixp,iyp,izp,ixp,iyp,izp+1,6)
c
           ttrans=tptmin
c 
c          try stencil B in 2D as a headwave
c
           if(max(time(ixp+1,iyp,izp+1),time(ixp,iyp+1,izp+1),
     +            time(ixp+1,iyp+1,izp+1)).lt.tptmin) 
     +       call timeb2d(ixp+1,iyp+1,izp+1,ixp,iyp+1,izp+1,
     +                  ixp+1,iyp,izp+1,ixp,iyp,izp+1,7)
           if(max(time(ixp-1,iyp,izp+1),time(ixp,iyp+1,izp+1),
     +            time(ixp-1,iyp+1,izp+1)).lt.tptmin) 
     +       call timeb2d(ixp-1,iyp+1,izp+1,ixp,iyp+1,izp+1,
     +                  ixp-1,iyp,izp+1,ixp,iyp,izp+1,7)
           if(max(time(ixp+1,iyp,izp+1),time(ixp,iyp-1,izp+1),
     +            time(ixp+1,iyp-1,izp+1)).lt.tptmin) 
     +       call timeb2d(ixp+1,iyp-1,izp+1,ixp,iyp-1,izp+1,
     +                  ixp+1,iyp,izp+1,ixp,iyp,izp+1,7)
           if(max(time(ixp-1,iyp,izp+1),time(ixp,iyp-1,izp+1),
     +            time(ixp-1,iyp-1,izp+1)).lt.tptmin) 
     +       call timeb2d(ixp-1,iyp-1,izp+1,ixp,iyp-1,izp+1,
     +                  ixp-1,iyp,izp+1,ixp,iyp,izp+1,7)
c
c          try stencil C in 1D as a headwave
c
           if(time(ixp+1,iyp,izp+1).lt.tptmin) 
     +       call time1d(ixp+1,iyp,izp+1,ixp,iyp,izp+1,8)
           if(time(ixp-1,iyp,izp+1).lt.tptmin) 
     +       call time1d(ixp-1,iyp,izp+1,ixp,iyp,izp+1,8)
           if(time(ixp,iyp+1,izp+1).lt.tptmin) 
     +       call time1d(ixp,iyp+1,izp+1,ixp,iyp,izp+1,8)
           if(time(ixp,iyp-1,izp+1).lt.tptmin) 
     +       call time1d(ixp,iyp-1,izp+1,ixp,iyp,izp+1,8)
c
        return      
c
      end
