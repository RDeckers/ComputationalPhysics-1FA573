set term term_type
set output output_file
set logscale xy
set key outside;
set key center top;
set format y "10^{%L}"
set format x "10^{%L}"
set ylabel "log absolute error"
set xlabel "h"
set grid;
plot\
 'data/table_1_1.dat' using 1:2 w lp title 'Symmetric 3-point',\
 'data/table_1_1.dat' using 1:3 w lp title 'Forward 2-point',\
 'data/table_1_1.dat' using 1:4 w lp title 'Backward 2-point'
set output