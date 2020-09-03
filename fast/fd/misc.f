c
c     version 1.1  Apr 1995
c
c     Misc routines for FD                
c
c        patched for g77 compilation
c
c        Drew Brenders (brenders@geoladm.geol.queensu.ca)
c
c     ----------------------------------------------------------------
c
      subroutine fnode(ix,iy,iz,torig,ttrans,ttmin,trmax,tthwmin,
     +              nreplace,treplace,hface,imodxz,imodxy,imodyz)
c
c     update things after timing one node
c
      include 'fd.par'
      include 'fd.com'
c
      integer hface(6)
c
      delt=torig-tptmin
      if(delt.gt.ttmin) then
        stencil(istencil)=stencil(istencil)+1
        if(ireverse.eq.1) then
          nreplace=nreplace+1
          treplace=treplace+delt
          if(delt.gt.trmax) trmax=delt
        end if
        if(istencil.gt.6) then
          if(ttrans-tptmin.gt.tthwmin)
     +      hwside(iside)=hwside(iside)+1
          if(ix.lt.hface(1)) hface(1)=ix
          if(ix.gt.hface(2)) hface(2)=ix
          if(iy.lt.hface(3)) hface(3)=iy
          if(iy.gt.hface(4)) hface(4)=iy
          if(iz.lt.hface(5)) hface(5)=iz
          if(iz.gt.hface(6)) hface(6)=iz
        end if
        ncalc=ncalc+1
        time(ix,iy,iz)=tptmin
        if(tptmin.lt.tminnew) tminnew=tptmin
c
        if(noplot.eq.0) then
          if(istencil.eq.3.or.istencil.gt.4) ipltnode=6
          if(istencil.gt.6) ipltnode=8
          if(inode.ne.1.or.ipltnode.ne.5) 
     +      call plotnode(ix,iy,iz,imodxz,imodxy,imodyz)
        end if
c
        if(iwrite.ne.1.or.ireverse.eq.1) return
c
        if(mod(ncalc,n10).eq.0) then
          write(io,5) nint(100.*float(ncalc)/float(nodestotal))
5         format(i5)
        end if
      end if
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine expand(iorder,istop,ioc,tmincurrent,nt,
     +                  nrepmem,nreplace,nfrev)
c
c     determine which face of the expanding cube to time next
c
      include 'fd.par'
      include 'fd.com'
c
      integer order(6)
      data order/5,6,3,4,1,2/
c
      if(iorder.eq.0) then
        tmincurrent=1.e11
        if(ireverse.eq.0.or.istop.ne.1) then
          do 2110 i=1,6
             if(side(i).ne.sidelimit(i).and.tminside(i).lt.
     +       tmincurrent) then
               tmincurrent=tminside(i)
               iside=i
             end if
2110      continue
        else
          do 3110 i=1,6
             if(side(i).ne.sidestop(i).and.tminside(i).lt.
     +       tmincurrent) then
               tmincurrent=tminside(i)
               iside=i
             end if
3110      continue
        end if
      else
        if(ireverse.eq.0.or.istop.ne.1) then
2120      ioc=ioc+1
          if(ioc.lt.7) then
            iside=order(ioc)
            if(side(iside).ne.sidelimit(iside)) go to 2130
          else
            ioc=0
          end if
          go to 2120
        else
2140      ioc=ioc+1
          if(ioc.lt.7) then
            iside=order(ioc)
            if(side(iside).ne.sidestop(iside)) go to 2130
          else
            ioc=0
          end if
          go to 2140
        end if
      end if
c
2130  continue
c
      tminnew=1.e10
      nt=0
      nrepmem=nreplace
      nfrev=nfrev+1
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine pass(nreplace,treplace,trmax,ncrev,reverse,nreverse,
     +                nhwmin,trmin,istop,hface,hpface,nfrev,irev)
c
c     after a single pass through the model, output headwaves and times
c     replaced and determine if another propagation is necessary
c
      include 'fd.par'
      include 'fd.com'
c
      integer reverse,revopp(6),hface(6),hpface(6)
      character face(6)*6
      data revopp/2,1,4,3,6,5/
     +     face/'  left',' right','  back',' front','   top','bottom'/
c
      nhw=hwside(1)+hwside(2)+hwside(3)+hwside(4)+hwside(5)+
     +    hwside(6)
      if(istop.eq.1) nhw=nhw*nodestotal/((side(2)-side(1)+1)*
     +               (side(4)-side(3)+1)*(side(6)-side(5)+1))
      if(iwrite.eq.1) then
        if(nhw.gt.0) then
          write(io,1005) (hwside(i),i=1,6)
1005      format(/'headwaves used:'/
     +         'left: ',i6,'   right:  ',i6/
     +         'back: ',i6,'   front:  ',i6/
     +         'top:  ',i6,'   bottom: ',i6)
        else
          write(io,1035)
1035      format(/'* no headwaves used *')
        end if
        if(ireverse.eq.1) then
          if(nreplace.gt.0) then
            write(io,1015) nreplace,treplace/nreplace,trmax
1015        format('number of traveltimes replaced: ',i7/
     +             'average time decrease:',f8.4,
     +          '   maximum time decrease: ',f8.4)
          else
            write(io,1025) nreplace
1025        format('number of traveltimes replaced: ',i7)
          end if
        end if
      end if
c
      if(ncrev.lt.nreverse.and.reverse.ne.0.and.nhw.gt.nhwmin.
     +  and.trmax.gt.trmin) then
        ncrev=ncrev+1
        if(ncrev.gt.1) then
          reverse=revopp(abs(reverse))
        else
          if(istop.eq.1) then
            do 1090 i=1,6
1090           sidestop(i)=side(i)
          end if
          ir6=max(6,abs(reverse))
          if((hwside(ir6).eq.0.and.reverse.gt.-7.
     +    and.reverse.lt.0).or.reverse.eq.7) then
            nhwmax=0
            do 1070 i=1,6
               if(hwside(i).gt.nhwmax) then
                 reverse=i
                 nhwmax=hwside(i)
               end if
1070        continue
          end if
        end if
        do 1020 i=1,6
1020       hwside(i)=0
        nside=5
        nstop=5
        reverse=abs(reverse)
        if(istop.ne.1) then
          side(1)=sidelimit(1)
          side(2)=sidelimit(2)
          side(3)=sidelimit(3)
          side(4)=sidelimit(4)
          side(5)=sidelimit(5)
          side(6)=sidelimit(6)
        else
          side(1)=sidestop(1)
          side(2)=sidestop(2)
          side(3)=sidestop(3)
          side(4)=sidestop(4)
          side(5)=sidestop(5)
          side(6)=sidestop(6)
        end if
        if(reverse.eq.1) then
          side(2)=hface(1)
          iface=hface(1)
        end if
        if(reverse.eq.2) then
          side(1)=hface(2)
          iface=hface(2)
        end if
        if(reverse.eq.3) then
          side(4)=hface(3)
          iface=hface(3)
        end if
        if(reverse.eq.4) then
          side(3)=hface(4)
          iface=hface(4)
        end if
        if(reverse.eq.5) then
          side(6)=hface(5)
          iface=hface(5)
        end if
        if(reverse.eq.6) then
          side(5)=hface(6)
          iface=hface(6)
        end if
        if(iwrite.eq.1) write(io,85) ncrev,face(reverse),iface
85      format(/'reverse propagation',i3,' from: ',a6/
     +,         'starting face: ',i5)
        hpface(1)=hface(1)
        hpface(2)=hface(2)
        hpface(3)=hface(3)
        hpface(4)=hface(4)
        hpface(5)=hface(5)
        hpface(6)=hface(6)
        hface(1)=nx
        hface(2)=1
        hface(3)=ny
        hface(4)=1
        hface(5)=nz
        hface(6)=1
        ireverse=1
        nreplace=0
        treplace=0.
        trmax=0.
        nfrev=0
        if(iclear.gt.0) icflag=iclear
        irev=1
      else
        irev=0
      end if
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine fside(nfrev,nreplace,nrepmem,incr,ip,ihpf,irep)
c
c     update things after timing one side of the cube
c
      include 'fd.par'
      include 'fd.com'
c
      irep=1
      side(iside)=side(iside)+incr
      if(side(iside).eq.sidelimit(iside)) nside=nside+1
      if(side(iside).eq.sidestop(iside)) nstop=nstop+1
      tminside(iside)=tminnew
c
      if(ireverse.eq.1) then
        if(nrepmem.eq.nreplace.and.ip.lt.ihpf) irep=0     
        if(iwrite.eq.1.and.irep.eq.0) write(io,785) side(iside)
785     format('     last face:     ',i5)
      end if
c
      return
      end
c
c     ----------------------------------------------------------------
c
      SUBROUTINE sort(arr,irr,n)
c
c     sort the elements of array arr in order of increasing size using
c     the quick sort technique
c
      INTEGER n,M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK),irr(n),b
      REAL a,temp
c
      do 300 i=1,n
         irr(i)=i
300    continue
c
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          b=irr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
            irr(i+1)=irr(i)
11        continue
          i=0
2         arr(i+1)=a
          irr(i+1)=b
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        itemp=irr(k)
        irr(k)=irr(l+1)
        irr(l+1)=itemp
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
          itemp=irr(l+1)
          irr(l+1)=irr(ir)
          irr(ir)=itemp
        endif
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
          itemp=irr(l)
          irr(l)=irr(ir)
          irr(ir)=itemp
        endif
        if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
          itemp=irr(l+1)
          irr(l+1)=irr(l)
          irr(l)=itemp
        endif
        i=l+1
        j=ir
        a=arr(l)
        b=irr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        itemp=irr(i)
        irr(i)=irr(j)
        irr(j)=itemp
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
        irr(l)=irr(j)
        irr(j)=b
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in sort2'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
