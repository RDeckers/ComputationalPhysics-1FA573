set term term_type
set output output_file
#unset key
# Axes
set ylabel "Θ"
set xlabel "b"

set format x '%.1f'
set format y '%.1f'
set format z '%.2f'

set multiplot;

set xrange [0:1];
set yrange [-3.16:3.16];

plot for [IDX=0:24:2] 'data/01_02.dat' i IDX using 1:2 w l lw 1.0 lc rgb 'grey' title "";
plot for [IDX=1:24:2] 'data/01_02.dat' i IDX using 1:2 w l lw .5 lc rgb 'black' title "";
plot for [IDX=0:24:6] 'data/01_02.dat' i IDX using 1:2 w l lw 1.9 title columnheader(1);


unset multiplot
set output
