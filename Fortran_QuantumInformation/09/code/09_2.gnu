# launch on terminal:
# gnuplot -e "ofile='outputfile.pdf'" -e "K=level" 09.gnu

# output file
set terminal pdfcairo enhanced font 'CMU Serif,18' size 6,4.5
set termopt enhanced
set encoding utf8
set output '../report/images/'.ofile

# set logscale for x and y axes
#set logscale x
#set logscale y

# format x,y axes ticks
set format x "%.1f"
set format y "%.1f"

# plot title
set title "Ising Haminltonian: {/:Italic k} = ".K." normalized eigenvalue"

# plot labels
set xlabel "Interaction strength λ"
set ylabel "Normalized {/:Italic E}_{/:Italic k}"

# set legend (key)
set key bottom left reverse Left box samplen 1

# fit curves
set fit quiet
set fit logfile '/dev/null'

# set colors palette inferno
set style line  1 lt 1 lc rgb '#000004' # black
set style line  2 lt 1 lc rgb '#1f0c48' # dark purple
set style line  3 lt 1 lc rgb '#550f6d' # dark purple
set style line  4 lt 1 lc rgb '#88226a' # purple
set style line  5 lt 1 lc rgb '#a83655' # red-magenta
set style line  6 lt 1 lc rgb '#e35933' # red
set style line  7 lt 1 lc rgb '#f9950a' # orange
set style line  8 lt 1 lc rgb '#f8c932' # yellow-orange
set style line  9 lt 1 lc rgb '#fcffa4' # light yellow

set grid
set grid xtics  lw 1
set grid ytics  lw 1
set grid mxtics lw 0.5
set grid mytics lw 0.5

# plot the graphic
plot for [C=2:9] "results/ising_".K."_K.dat" u 1:C with lp ls (C-1) pt 7 ps 0.35 lw 2 title "{/:Italic N} = ".(C+1)."  "
