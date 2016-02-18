set term term_type
set output output_file
set logscale xy
set key outside;
set key center top;
set format y "10^{%L}"
set format x "10^{%L}"
set ylabel "solution at r=0"
set xlabel "h"
set grid;
plot\
 'data/exercise_3_2.dat' using 1:2 w lp title 'Numerov'
set output
