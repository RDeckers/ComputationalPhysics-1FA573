set term term_type
set output output_file
#unset key
# Axes
set ylabel "|error(Θ)|"
set xlabel "b"

set format x '%.1f'
set format y "10^{%L}"
set logscale y

plot 'data/01_03.dat' using 1:3 w l title "Partially approximated",\
  '' u 1:4 w l title "Fully approximated",\
  '' u 1:5 w l title "Realistically approximated";
set output
