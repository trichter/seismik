#!/bin/tcsh

set imax=`tail -1 i.dat | awk '{print $1}'`

grd2xyz R$imax.grd > dummy1
grd2xyz V$imax.grd > dummy2

paste dummy2 dummy1 | awk '{print $1, 285-$2, $3, $6}' > dummy3


cat dummy3 | awk '{if($3<5.5)print $0}'  > VR$imax.dat


#\rm dummy1 dummy2 dummy3

exit
