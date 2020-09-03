#!/bin/csh -f
#
# This will make fd, nray, xray, and all other FAST programs
#
# Drew Brenders (brenders@geoladm.geol.queensu.ca)
#
#
echo ""
echo "Compiling FAST!  Note compilation warnings on stdout!"
echo ""
echo "***************************************************"

cd fd
make
cd ../ray
make nray
make xray
cd ../zslice
make
cd ../misc
cp ../useful_tools/compile_misc.csh .
./compile_misc.csh
cd ..

echo ""
echo "***************************************************"
echo ""
echo "Removing objects..."
echo ""
rm -Rf ./fd/*.o ./ray/*.o ./zslice/*.o 
echo "***************************************************"
echo ""
echo "Finished!  Compiled binaries are in ./bin"
echo ""
