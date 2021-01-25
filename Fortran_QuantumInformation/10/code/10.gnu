# launch on terminal:
# gnuplot -e "ofile='outputfile.pdf'" -e "N_min=N_min" -e "N_max=N_max" -e "I=iterations" 10.gnu

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
set title "RSRG algorithm results for ground state"

# plot labels
set xlabel "Interaction strength λ"
set ylabel "{/:Italic E}_{gs}(λ)"

# set legend (key)
set key bottom center reverse Left box samplen 1

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
set yrange [-4:0]

mf(x) = abs(x)>2 ? -abs(x) : -1.0 -x**2/4.0

# plot the graphic
plot for [k=N_min:N_max] "results/ising_RSRG_".I."_I.dat" u 1:k with lp ls (k-N_min+1)     pt 7 ps 0.35 lw 2 title "{/:Italic N} = ".k, \
                         mf(x)                                  with l  ls (N_max-N_min+2) dashtype 4   lw 2 title "Mean field solution"
