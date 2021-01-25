# launch on terminal:
# gnuplot -e "ofile='outputfile.pdf'" 05d.gnu

# output file
set terminal pdfcairo font 'CMU Serif,18' size 6,4.5
set termopt enhanced
set output '../report/images/'.ofile

# format x,y axes ticks
set format x "%.1f"
set format y "%.1f"

# plot title
set title "Normalized spacing between eigenvalues"

# plot labels
set xlabel "s"
set ylabel "P(s)"

# set legend (key)
set key top right reverse Left box samplen 1

# fit curves
set fit quiet
set fit logfile '/dev/null'

f1(x) = a1*(x**b1)*exp(-c1*x**d1)
f2(x) = a2*(x**b2)*exp(-c2*x**d2)
f3(x) = a3*(x**b3)*exp(-c3*x**d3)
f4(x) = a4*(x**b4)*exp(-c4*x**d4)

fit f1(x) "results/herm_50_samples_2500_N_glob.dat"          u 1:2 via a1, b1, c1, d1
fit f2(x) "results/herm_50_samples_2500_N_250_level_loc.dat" u 1:2 via a2, b2, c2, d2
fit f3(x) "results/herm_50_samples_2500_N_50_level_loc.dat"  u 1:2 via a3, b3, c3, d3
fit f4(x) "results/herm_50_samples_2500_N_10_level_loc.dat"  u 1:2 via a4, b4, c4, d4

set fit errorvariables

# print fit results on stdout
# global
print sprintf('Fit results (global):')
print sprintf('a1=%2.3f \pm %2.3f', a1, a1_err)
print sprintf('b1=%2.3f \pm %2.3f', b1, b1_err)
print sprintf('c1=%2.3f \pm %2.3f', c1, c1_err)
print sprintf('d1=%2.3f \pm %2.3f', d1, d1_err)
# local: level 100
print sprintf('Fit results (local, level=250):')
print sprintf('a2=%2.3f \pm %2.3f', a2, a2_err)
print sprintf('b2=%2.3f \pm %2.3f', b2, b2_err)
print sprintf('c2=%2.3f \pm %2.3f', c2, c2_err)
print sprintf('d2=%2.3f \pm %2.3f', d2, d2_err)
# local: level 50
print sprintf('Fit results (local, level=50):')
print sprintf('a3=%2.3f \pm %2.3f', a3, a3_err)
print sprintf('b3=%2.3f \pm %2.3f', b3, b3_err)
print sprintf('c3=%2.3f \pm %2.3f', c3, c3_err)
print sprintf('d3=%2.3f \pm %2.3f', d3, d3_err)
# local: level 10
print sprintf('Fit results (local, level=10):')
print sprintf('a4=%2.3f \pm %2.3f', a4, a4_err)
print sprintf('b4=%2.3f \pm %2.3f', b4, b4_err)
print sprintf('c4=%2.3f \pm %2.3f', c4, c4_err)
print sprintf('d4=%2.3f \pm %2.3f', d4, d4_err)

# set grid
set grid xtics  lw 1
set grid ytics  lw 1
set grid mxtics lw 0.5
set grid mytics lw 0.5

#plot the graphic
plot "results/herm_50_samples_2500_N_glob.dat"          u 1:2 with p ls 1 lc rgb '#1f77b4'      pt 7 lw 0 ps 0.5 title "P(s) (global)", \
     f1(x)                                                    with l ls 1 lc rgb '#1f77b4' lt 1      lw 2        title sprintf("P(s)=%.2f s^{%.2f} e^{-%.2f s^{%.2f}}", a1, b1, c1, d1), \
     "results/herm_50_samples_2500_N_250_level_loc.dat" u 1:2 with p ls 1 lc rgb '#ff7f0e'      pt 7 lw 0 ps 0.5 title "P(s) (local, level=250)", \
     f2(x)                                                    with l ls 1 lc rgb '#ff7f0e' lt 1      lw 2        title sprintf("P(s)=%.2f s^{%.2f} e^{-%.2f s^{%.2f}}", a2, b2, c2, d2), \
     "results/herm_50_samples_2500_N_50_level_loc.dat"  u 1:2 with p ls 1 lc rgb '#2ca02c'      pt 7 lw 0 ps 0.5 title "P(s) (local, level=50)", \
     f3(x)                                                    with l ls 1 lc rgb '#2ca02c' lt 1      lw 2        title sprintf("P(s)=%.2f s^{%.2f} e^{-%.2f s^{%.2f}}", a3, b3, c3, d3), \
     "results/herm_50_samples_2500_N_10_level_loc.dat"  u 1:2 with p ls 1 lc rgb '#d62728'      pt 7 lw 0 ps 0.5 title "P(s) (local, level=10)", \
     f4(x)                                                    with l ls 1 lc rgb '#d62728' lt 1      lw 2        title sprintf("P(s)=%.2f s^{%.2f} e^{-%.2f s^{%.2f}}", a4, b4, c4, d4)
