# launch on terminal:
# gnuplot -e "ofile='outputfile.pdf'" 08_ns.gnu

# output file
set terminal pdfcairo enhanced font 'CMU Serif,18' size 6,4.5
set termopt enhanced
set encoding utf8
set output '../report/images/'.ofile

# set logscale for x and y axes
#set logscale x
set logscale y

# format x,y axes ticks
set format x  "%.0f"
set format y "10^{%L}"

# plot title
set title "CPU time needed for non-separable state initialization"

# plot labels
set xlabel "Number of subsystems N"
set ylabel "CPU time [s]"

# set legend (key)
set key top left reverse Left box samplen 1

# fit curves
set fit quiet
set fit logfile '/dev/null'

set style line  1 lt 1 lc rgb '#1f77b4' # blue
set style line  2 lt 1 lc rgb '#ff7f0e' # orange
set style line  3 lt 1 lc rgb '#2ca02c' # green
set style line  4 lt 1 lc rgb '#d62728' # red
set style line  5 lt 1 lc rgb '#9467bd'
set style line  6 lt 1 lc rgb '#8c564b'
set style line  7 lt 1 lc rgb '#e377c2'
set style line  8 lt 1 lc rgb '#7f7f7f'
set style line  9 lt 1 lc rgb '#bcbd22'
set style line 10 lt 1 lc rgb '#17becf'

set grid
set grid xtics  lw 1
set grid ytics  lw 1
set grid mxtics lw 0.5
set grid mytics lw 0.5

set xr [2:28]
# plot the graphic
plot "results/nonseparable_2_D_2.dat" u 1:2 with lp ls 1 pt 7 ps 0.35 lw 2 title "CPU time: d=2  "