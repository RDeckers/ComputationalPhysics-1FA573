set term term_type
set output output_file
set logscale y
set key outside;
set key center top;
set format y "10^{%L}"
set ylabel "log absolute error"
set xlabel "iteration"
set grid;
plot\
 'data/table_1_4.dat' using 1:2 w lp title 'search',\
 'data/table_1_4.dat' using 1:3 w lp title 'Newton-Rhapson',\
 'data/table_1_4.dat' using 1:4 w lp title 'secant'
unset output
