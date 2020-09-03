      program main     
      parameter (n2max=1700,n1max=700)
      integer i,j,n1,n2
      real sep(n2max,n1max)
      integer*2 fast(n2max,n1max)
      character*80 sepfile,fastfile     

      write(*,*)'####################################'
      write(*,*)'####### Bin sep to bin Fast ########'
      write(*,*)'####################################'
      write(*,'(/,"Enter SEP binary file name to reformat: ",$)')
      read(*,100) sepfile
      write(*,'(/,"Enter Fast binary file name to reformat: ",$)')
      read(*,100) fastfile
100    format(a80)
     

      write(*,'(/,"Enter nb of traces & nb of samples per trace",
     & /," N2 and N1 : ",$)')
      read(*,*) n2,n1

      open(10,file=fastfile,form='unformatted')
       open(11,file=sepfile,
     +     access='direct',form='unformatted',recl=4*n1)

      if (n1.gt.n1max) then
      	print *,'change n1max parameter an recompile!!'
      	stop
      endif
      if (n2.gt.n2max) then
      	print *,'change n2max parameter an recompile!!'
      	stop
      endif
      
      if(n1*n2.gt.1E+4) write(*,*)' ... I''m runing ...'

      
      do 3 i=1,n2
         read(11,rec=i) (sep(i,j),j=1,n1)
3     continue
         close(11)
         
      do j=1,n1  
        do i=1,n2 
        fast(i,j)=nint(sep(i,j)) 
        enddo
      enddo        
      do j=1,n1
      write (10) (fast(i,j),i=1,n2)
      enddo

        close(10)
         end

