# launch on terminal:
# gnuplot -e "k=k_eigf" -e "ofile='outputfile.pdf'" 06_ef.gnu
# gnuplot -e "k=1" -e "ofile='1_eigenfunction.pdf'" 06_ef.gnu
# gnuplot -e "k=2" -e "ofile='2_eigenfunction.pdf'" 06_ef.gnu
# gnuplot -e "k=3" -e "ofile='3_eigenfunction.pdf'" 06_ef.gnu
# gnuplot -e "k=4" -e "ofile='4_eigenfunction.pdf'" 06_ef.gnu

# output file
set terminal pdfcairo enhanced font 'CMU Serif,18' size 6,4.5
set termopt enhanced
set encoding utf8
set output '../report/images/'.ofile

# format x,y axes ticks
set format x "%.1f"
set format y "%.2f"

# plot title
set title "1D Harmonic Oscillator Eigenfunctions"

# plot labels
set xlabel "{/:Italic x}"
set ylabel "Ψ({/:Italic x})"
#set ylabel "{/Symbol y}({/:Italic x})"

# set legend (key)
set key top right reverse Left box samplen 1

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
plot [-5:5] "results/".k."_eigenfunction.dat" u 1:2 with l ls k lw 2 title "k=".k." eigenfunction"
