set term term_type
set output output_file
set logscale xy
set key outside;
set key center top;
set format y "10^{%L}"
set format x "10^{%L}"
set ylabel "absolute rel. error"
set xlabel "h"
set grid;
plot\
 'data/exercise_2_1.dat' using 1:2 w lp title 'Euler',\
 'data/exercise_2_1.dat' using 1:3 w lp title 'Taylor',\
 'data/exercise_2_1.dat' using 1:4 w lp title 'Implicit'
set output
