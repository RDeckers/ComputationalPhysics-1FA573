set term term_type
set output output_file
set logscale xy
set key outside;
set key center top;
set format y "10^{%L}"
set format x "10^{%L}"
set ylabel "absolute drift after 100 iterations"
set xlabel "dt"
set grid;
plot\
 'data/exercise_2_5.dat' using 1:2 w lp title 'p',\
 'data/exercise_2_5.dat' using 1:3 w lp title 'y',\
 'data/exercise_2_5.dat' using 1:(sqrt($2**2+$3**2)) w lp title 'sqrt(p^2+y^2)'
set output
