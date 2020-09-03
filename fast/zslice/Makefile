#
# makefile for ZSLICE
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
XZSLICE_OBJS= main.o misc.o blkdat.o $(PLTLIB)/pltsub.o $(PLTLIB)/xpltlib.o $(PLTLIB)/xbuplot.o
#
# Executable name
#
XZSLICE_EXEC=$(EXEDIR)/xzslice
#
xzslice: $(XZSLICE_EXEC)
$(XZSLICE_EXEC): ${XZSLICE_OBJS}
	$(FC) -o main ${XZSLICE_OBJS} -L/usr/X11R6/lib -lX11
	mv main ${XZSLICE_EXEC}
#
# dependencies
#
main.o     :  zslice.par zslice.com
misc.o     :  zslice.par zslice.com
blkdat.o   :  zslice.par zslice.com

$(PLTLIB)/pltsub.o: FORCE
	cd $(@D); $(MAKE) $(@F)
FORCE:

$(PLTLIB)/xpltlib.o: FORCE
	cd $(@D); $(MAKE) $(@F)
FORCE:

$(PLTLIB)/nopltlib.o: FORCE
	cd $(@D); $(MAKE) $(@F)
FORCE:

$(PLTLIB)/xbuplot.o: FORCE
	cd $(@D);cc -O -I/usr/X11R6/include/X11 -c xbuplot.c
FORCE:


