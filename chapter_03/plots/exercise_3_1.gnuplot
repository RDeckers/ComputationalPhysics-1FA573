set term term_type
set output output_file
set logscale xy
set key outside;
set key center top;
set format y "10^{%L}"
set format x "10^{%L}"
set ylabel "absolute error"
set xlabel "h"
set grid;
plot\
 'data/exercise_3_1.dat' using 1:2 w lp title 'Numerov',\
 'data/exercise_3_1.dat' using 1:3 w lp title 'RK4'
set output
