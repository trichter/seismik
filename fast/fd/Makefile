#
# makefile for FD
#
# patched for g77 compilation
#
# Drew Brenders (brenders@geoladm.geol.queensu.ca)
#
EXEDIR=../bin
PLTLIB=../pltlib
#
# compile flags
#
FC=f77 -static
FFLAGS =  -O3
#CFLAGS =  -O -I/usr/openwin/include
CFLAGS =  -O 
#
#rules for compiling .f files:
#
.f.o:
	$(COMPILE.f) $(FFLAGS) $< $(INCS)
.c.o:
	$(COMPILE.c) $(CFLAGS) $< $(INCS)
#---------------------------------------------------------------------------
#
# Source objects
#
NFD_OBJS= main.o model.o time.o findiff.o findiff2d.o stencils.o stencils2d.o misc.o plt.o blkdat.o $(PLTLIB)/nopltlib.o
#
# Executable name
#
NFD_EXEC=$(EXEDIR)/nfd
#
nfd: $(NFD_EXEC)
$(NFD_EXEC): ${NFD_OBJS}
	$(FC) -o main ${NFD_OBJS}
	mv main ${NFD_EXEC}
#
# dependencies
#
main.o       :  fd.par fd.com
model.o      :  fd.par fd.com
time.o       :  fd.par fd.com
findiff.o    :  fd.par fd.com
findiff2d.o  :  fd.par fd.com
stencils.o   :  fd.par fd.com
stencils2d.o :  fd.par fd.com
misc.o       :  fd.par fd.com
plt.o        :  fd.par fd.com
blkdat.o     :  fd.par fd.com
 
$(PLTLIB)/nopltlib.o: FORCE
	cd $(@D); $(MAKE) $(@F)
FORCE:

