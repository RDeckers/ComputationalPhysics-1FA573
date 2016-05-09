set term term_type
set output output_file
#unset key
# Axes
set ylabel "dσ/dΩ"
set xlabel "Θ"

#set format x '%.1f'
set logscale y
set xrange [0:3.14]
plot for [IDX=0:6] 'data/01_04b.dat' i IDX using 1:2 w l lw 1.6 title columnheader(1);
set output
