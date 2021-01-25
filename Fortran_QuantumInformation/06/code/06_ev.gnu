# launch on terminal:
# gnuplot -e "ofile='outputfile.pdf'"  06_ev.gnu
# gnuplot -e "ofile='eigenvalues.pdf'" 06_ev.gnu

# output file
set terminal pdfcairo enhanced font 'CMU Serif,18' size 6,4.5
set termopt enhanced
set encoding utf8
set output '../report/images/'.ofile

# format x,y axes ticks
# set format x "%.0f"
# set format y "%.2f"
set format x "10^{%L}"
set format y "10^{%L}"

# plot title
# set title "1D Harmonic Oscillator Eigenvalues: |{/:Italic E}_{{/:Italic k},num} - {/:Italic E}_{{/:Italic k},th}| / {/:Italic E}_{{/:Italic k},th}"
set title "1D Harmonic Oscillator Eigenvalues: RelErr({/:Italic E}_{{/:Italic k},num},{/:Italic E}_{{/:Italic k},th}) [%]"

# plot labels
set xlabel "{/:Italic k}"
# set ylabel "|{/:Italic E}_{{/:Italic k},num} - {/:Italic E}_{{/:Italic k},th}| / {/:Italic E}_{{/:Italic k},th}"
set ylabel "RelErr({/:Italic E}_{{/:Italic k},num},{/:Italic E}_{{/:Italic k},th}) [%]"

# set legend (key)
set key top left reverse Left box samplen 1

# set logscale for x and y axes
set logscale x
set logscale y

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

# set grid
set grid xtics  lw 1
set grid ytics  lw 1
set grid mxtics lw 0.5
set grid mytics lw 0.5

# plot the graphic
# plot for [i=1:5] "results/".i."_eigenfunction.dat" u 1:2 with l ls i lw 2 title "k=".i." eigenfunction"
plot "results/eigenvalues.dat" u 1:(100*abs($4)/column(2)) with lp ls 1 pt 7 ps 0.35 lw 2 title "RelErr({/:Italic E}_{{/:Italic k},num},{/:Italic E}_{{/:Italic k},th}) [%]"
