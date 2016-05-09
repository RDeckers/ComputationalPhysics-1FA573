set term term_type
set output output_file
#unset key
# Axes
set ylabel "Θ"
set xlabel "b"

set format x '%.1f'
set format y "%.1f"

plot for [IDX=0:5] 'data/01_04.dat' i IDX using 1:2 w l title columnheader(1);

set output
