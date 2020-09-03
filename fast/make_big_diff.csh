#!/bin/csh -f

set ori_dir = ../fast+brender

rm fast_g77_big+marco.diff
set list = `ls fd/Make* fd/*.f fd/*.par fd/*.doc fd/*.c fd/*.com`
foreach file ($list)
 echo "diff -u -r $ori_dir/$file $file > fast_g77_big+marco.diff" >> fast_g77_big+marco.diff
 diff -u -r $ori_dir/$file $file >> fast_g77_big+marco.diff
end
set list = `ls misc/Make* misc/*.f misc/*.par misc/*.doc misc/*.c misc/*.com`
foreach file ($list)
 echo "diff -u -r $ori_dir/$file $file > fast_g77_big+marco.diff" >> fast_g77_big+marco.diff
 diff -u -r $ori_dir/$file $file >> fast_g77_big+marco.diff
end
set list = `ls pltlib/Make* pltlib/*.f pltlib/*.par pltlib/*.doc pltlib/*.c pltlib/*.com`
foreach file ($list)
 echo "diff -u -r $ori_dir/$file $file > fast_g77_big+marco.diff" >> fast_g77_big+marco.diff
 diff -u -r $ori_dir/$file $file >> fast_g77_big+marco.diff
end
set list = `ls ray/Make* ray/*.f ray/*.par ray/*.doc ray/*.c ray/*.com`
foreach file ($list)
 echo "diff -u -r $ori_dir/$file $file > fast_g77_big+marco.diff" >> fast_g77_big+marco.diff
 diff -u -r $ori_dir/$file $file >> fast_g77_big+marco.diff
end
set list = `ls zslice/Make* zslice/*.f zslice/*.par zslice/*.doc zslice/*.c zslice/*.com`
foreach file ($list)
 echo "diff -u -r $ori_dir/$file $file > fast_g77_big+marco.diff" >> fast_g77_big+marco.diff
 diff -u -r $ori_dir/$file $file >> fast_g77_big+marco.diff
end
