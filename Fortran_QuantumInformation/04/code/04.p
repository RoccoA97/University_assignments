# launch on terminal:
# gnuplot -e "titlename='titlename'" -e "mean_col=mean_col_number" -e "std_col=std_col_number" 04.p

# output file
set terminal pdfcairo font 'CMU Serif,11' size 4,3
set output '../report/images/flag'.titlename.'.pdf'

# set logscale for x and y axes
set logscale x
set logscale y

# format y axis ticks
set format y "10^{%L}"

# plot title
set title titlename." optimization"#'Performance (-O0 flag)'

# plot labels
set xlabel "Matrix size n"
set ylabel "CPU time [s]"

# Set linestyle 1 to blue (#0060ad)
set colorsequence podo
set style line 1 \
    linecolor rgb '#0060ad' \
    linetype 1 linewidth 2 \
    pointtype 7 pointsize 0.5

set key top left reverse Left box samplen 1

# fit curves
set fit quiet
set fit logfile '/dev/null'

f(x) = a1*x + b1
g(x) = a2*x + b2
h(x) = a3*x + b3

#[log(200):]
#[log(200):]
#[log(200):]

fit f(x) 'matmul_col.dat' u (log($1)):(log(column(mean_col))):(log(column(mean_col))) yerrors via a1,b1
fit g(x) 'matmul_row.dat' u (log($1)):(log(column(mean_col))):(log(column(mean_col))) yerrors via a2,b2
fit h(x) 'matmul.dat'     u (log($1)):(log(column(mean_col))):(log(column(mean_col))) yerrors via a3,b3

set fit errorvariables

print sprintf('Fit result:  a1=%2.3f \pm %2.3f', a1, a1_err)
print sprintf('Fit result:  b1=%2.3f \pm %2.3f', b1, b1_err)
print sprintf('Fit result:  a2=%2.3f \pm %2.3f', a2, a2_err)
print sprintf('Fit result:  b2=%2.3f \pm %2.3f', b2, b2_err)
print sprintf('Fit result:  a3=%2.3f \pm %2.3f', a3, a3_err)
print sprintf('Fit result:  b3=%2.3f \pm %2.3f', b3, b3_err)

# set xrange
# set xrange [40:2500]

set grid xtics lw 1
set grid ytics lw 1
set grid mxtics lw 0.5
set grid mytics lw 0.5

#plot the graphic
plot 'matmul_col.dat' u 1:mean_col:std_col with yerrorbars  ls 1 lc 1 lw 0 title "matmul\\_col", \
     'matmul_row.dat' u 1:mean_col:std_col with yerrorbars  ls 1 lc 2 lw 0 title "matmul\\_row", \
     'matmul.dat'     u 1:mean_col:std_col with yerrorbars  ls 1 lc 3 lw 0 title "matmul", \
     exp(f(log(x)))                        with linespoints ls 1 lc 1 lw 2 pointsize 0 title sprintf("log(y)=%.2f log(x) %+.2f", a1, b1), \
     exp(g(log(x)))                        with linespoints ls 1 lc 2 lw 2 pointsize 0 title sprintf("log(y)=%.2f log(x) %+.2f", a2, b2), \
     exp(h(log(x)))                        with linespoints ls 1 lc 3 lw 2 pointsize 0 title sprintf("log(y)=%.2f log(x) %+.2f", a3, b3)
     #"log(y) = ${a3}*log(x) + ${b3}"#matmul fit"
