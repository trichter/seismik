#! /bin/csh -f


set COMP="gfortran"
set FLAG_MAIN=" -fbounds-check -mcmodel=medium"
set FLAG_SUB=" -c -w -fbounds-check  -mcmodel=medium "
set bin="bin2d"

if ( ! -d $bin )then
mkdir $bin
endif

\rm $bin/*

cd pltlib
set list=`ls *.f`
foreach file (${list})
  echo ${COMP} ${file} ${FLAG_SUB}
  ${COMP} ${file} ${FLAG_SUB}
end
gcc -c xbuplot.c
cd ../


cd fd
set list=`ls *.f`
set list_obj=`echo $list | sed 's/\.f/\.o/g'`
foreach file (${list})
  echo ${COMP} ${file} ${FLAG_SUB}
  ${COMP} ${file} ${FLAG_SUB}
end
${COMP} -o nfd  ${list_obj} ../pltlib/nopltlib.o ${FLAG_MAIN}
mv nfd ../${bin}
cd ../

cd ray
set list=`ls *.f`
set list_nray="main.o time.o model.o tomo.o plt.o blkdat.o ../pltlib/nopltlib.o"
set list_xray="main.o time.o model.o tomo.o plt.o blkdat.o ../pltlib/pltsub.o ../pltlib/xpltlib.o ../pltlib/xbuplot.o"

foreach file (${list})
  echo ${COMP} ${file} ${FLAG_SUB}
  ${COMP} ${file} ${FLAG_SUB}
end
${COMP} -o nray  ${list_nray} ${FLAG_MAIN}
mv nray ../${bin}
${COMP} -o xray  ${list_xray} ${FLAG_MAIN}
mv xray ../${bin}
cd ../


cd misc
set list=`ls *.f`
foreach file (${list})
  echo ${COMP} ${file} -o ${file:r} ${FLAG_MAIN}
  ${COMP} ${file} -o ${file:r} ${FLAG_MAIN}
  mv ${file:r} ../${bin}
end
cd ../



cd misc
set list=`ls *.f`
set list_obj=`echo $list | sed 's/\.f/\.o/g'`
foreach file (${list})
  echo ${COMP} ${file} ${FLAG_SUB}
  ${COMP} ${file} ${FLAG_SUB}
end
foreach file (${list_obj})
  echo ${COMP} ${file} -o ${file:r} 
  ${COMP} ${file} -o ${file:r}
  mv ${file:r} ../${bin}
end
cd ../


exit
foreach file (${list})
  echo ${COMP} ${file} -o ${file:r} 
  ${COMP} ${file} -o ${file:r}
  mv ${file:r} ../${bin}
end

cd ../
