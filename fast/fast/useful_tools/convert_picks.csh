#!/bin/csh -f
set list = `ls ./tmp/fd*`
foreach file ($list)
echo $file > convert.in
echo ${file:r}.picks >> convert.in
~/fast_g77/bin/rec_binary < convert.in
rm -f convert.in
end
