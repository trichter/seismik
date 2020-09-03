#!/bin/csh -f
set list = `ls vel*.asc`
foreach file ($list)
echo 25 > convert.in
echo 1000 >> convert.in
transp < $file > ${file}.transp n1=601 nbpe=6
ascii2su ${file}.transp < convert.in
su2bin ${file}.transp.su
echo ${file}.transp.su > convert.in
echo ${file:r} >> convert.in
echo "601 25" >> convert.in
sep2fast < convert.in
rm -f convert.in
rm -f ${file}.transp
rm -f ${file}.transp.su
rm -f ${file}.transp.su.bin
end
