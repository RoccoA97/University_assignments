# launch on terminal:
# gnuplot -e "ofile='outputfile.pdf'" 07_expv.gnu

# output file
set terminal pdfcairo enhanced font 'CMU Serif,18' size 6,4.5
set termopt enhanced
set encoding utf8
set output '../report/images/'.ofile

# format x,y axes ticks
set format x  "%.1f"
set format y  "%.1f"
set format cb "%.1f"

# plot title
set title "Position expectation value over time"

# plot labels
set xlabel "{/:Italic x} [a.u.]"
set ylabel "{/:Italic t} [a.u.]"
#set ylabel "{/Symbol y}({/:Italic x})"

# set legend (key)
set key top left reverse Left box samplen 1

# fit curves
set fit quiet
set fit logfile '/dev/null'


# set style line  1 lt 1 lc rgb '#440154' # dark purple
# set style line  2 lt 1 lc rgb '#472c7a' # purple
# set style line  3 lt 1 lc rgb '#3b518b' # blue
# set style line  4 lt 1 lc rgb '#2c718e' # blue
# set style line  5 lt 1 lc rgb '#21908d' # blue-green
# set style line  6 lt 1 lc rgb '#27ad81' # green
# set style line  7 lt 1 lc rgb '#5cc863' # green
# set style line  8 lt 1 lc rgb '#aadc32' # lime green
# set style line  9 lt 1 lc rgb '#fde725' # yellow

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

# plot the graphic
# plot for [i=1:5] "results/".i."_eigenfunction.dat" u 1:2 with l ls i lw 2 title "k=".i." eigenfunction"
# plot [-5:5] "results/".k."_eigenfunction.dat" u 1:2 with l ls k lw 2 title "k=".k." eigenfunction"
# set xrange [-3:3]
# set yrange [0:2]
stats  "results/wvfc_cmap.dat"  u 1 nooutput
x_min = STATS_min
x_max = STATS_max
stats  "results/wvfc_cmap.dat"  u 2 nooutput
y_min = STATS_min
y_max = STATS_max

# set xr [STATS_min_x:STATS_max_x]
# set yr [STATS_min_y:STATS_max_y]
set xr [x_min:x_max]
set yr [y_min:y_max]

plot "results/wvfc_expv.dat" u 2:1 with l ls 1      lw 4 title "⟨{/:Italic x}⟩", \
     "results/wvfc_expv.dat" u 3:1 with l ls 2 dt 4 lw 3 title "⟨{/:Italic x}⟩ ± ⟨Δ{/:Italic x}⟩", \
     "results/wvfc_expv.dat" u 4:1 with l ls 2 dt 4 lw 3 notitle
