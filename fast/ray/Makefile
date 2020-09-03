#
# makefile for RAY
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
XRAY_OBJS= main.o time.o model.o tomo.o plt.o blkdat.o $(PLTLIB)/pltsub.o $(PLTLIB)/xpltlib.o $(PLTLIB)/xbuplot.o
NRAY_OBJS= main.o time.o model.o tomo.o plt.o blkdat.o $(PLTLIB)/nopltlib.o
#
# Executable name
#
XRAY_EXEC=$(EXEDIR)/xray
NRAY_EXEC=$(EXEDIR)/nray
#
xray: $(XRAY_EXEC)
$(XRAY_EXEC): ${XRAY_OBJS}
	$(FC) -o main ${XRAY_OBJS} -L/usr/X11R6/lib -lX11
	mv main ${XRAY_EXEC}
#
nray: $(NRAY_EXEC)
$(NRAY_EXEC): ${NRAY_OBJS}
	$(FC) -o main ${NRAY_OBJS} 
	mv main ${NRAY_EXEC}
#
# dependencies
#
main.o     :  ray.par ray.com
time.o     :  ray.par ray.com
model.o    :  ray.par ray.com
tomo.o     :  ray.par ray.com
plt.o      :  ray.par ray.com
blkdat.o   :  ray.par ray.com
 
$(PLTLIB)/pltsub.o: FORCE
	cd $(@D); $(MAKE) $(@F)
FORCE:
 
$(PLTLIB)/xpltlib.o: FORCE
	cd $(@D); $(MAKE) $(@F)
FORCE:
 
$(PLTLIB)/xbuplot.o: FORCE
	cd $(@D);cc -O -I/usr/X11R6/include/X11 -c xbuplot.c
FORCE:

$(PLTLIB)/nopltlib.o: FORCE
	cd $(@D); $(MAKE) $(@F)
FORCE:

