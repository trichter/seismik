#! /bin/csh -f

set list=`ls *.f`
set FLAG="-o"
echo $list
foreach file (${list})
  echo $file
  f77 ${file} ${FLAG} ${file:r}
  mv ${file:r} ../bin/.

end
