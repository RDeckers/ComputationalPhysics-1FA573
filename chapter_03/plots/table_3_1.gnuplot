set term term_type
set output output_file
set key outside;
set key center top;
set logscale y;
set format y "10^{%L}"
set ylabel "absolute rel. error"
set xlabel "r"
set grid;
plot\
 'data/table_3_1.dat' using 1:2 w lp title 'exact y1',\
 'data/table_3_1.dat' using 1:3 w lp title 'perturbed y1',\
 'data/table_3_1.dat' using 1:4 w lp title 'linear corrected'
set output
