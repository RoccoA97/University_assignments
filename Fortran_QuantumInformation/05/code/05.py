import os
import subprocess
import math





def run_fortran(RUN, type, samples, N, nbins, xmin, xmax, g_or_l, level):
    print("Doing:")
    print("\ttype    = " + type)
    print("\tsamples = " + str(samples))
    print("\tN       = " + str(N))
    print("\tnbins   = " + str(nbins))
    print("\txmin    = " + str(xmin))
    print("\txmax    = " + str(xmax))
    print("\tlocal   = " + str(g_or_l))
    print("\tlevel   = " + str(level))
    if RUN:
        subprocess.run([exec, type, str(samples), str(N), str(nbins), str(xmin), str(xmax), str(g_or_l), str(level)])
    print("\n")


def run_gnuplot(RUN, type, samples, N, gnufile):
    print("\nNow fitting and plotting:")
    cmd  = "gnuplot -e "
    cmd += "\"ofile='" + type +  "_" + str(samples) + "_samples_" + str(N) + "_N.pdf" + "'\" "
    cmd += gnufile
    if RUN:
        os.system(cmd)
    return



RUN_FORT = False
RUN_GNU  = True
comp     = "gfortran"
c_flags  = "-Wall -ffree-line-length-0"
exec     = "./05.o"
src      = "05.f90"
lapack   = "-L/usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0 -llapack"

samples = 50
herm    = "herm"
diag    = "diag"
N       = 2500
nbins   = 100
xmin    = 0
xmax    = 5
loc     = 1
glob    = 0

# compile src
os.system(' '.join([comp, c_flags, "-o", exec, src, lapack]))

# run global and localversions
# global
run_fortran(RUN=RUN_FORT, type="herm", samples=samples, N=N, nbins=nbins, xmin=xmin, xmax=xmax, g_or_l=glob, level=N)
run_fortran(RUN=RUN_FORT, type="diag", samples=samples, N=N, nbins=nbins, xmin=xmin, xmax=xmax, g_or_l=glob, level=N)
# local: levels (10,50,100)
run_fortran(RUN=RUN_FORT, type="herm", samples=samples, N=N, nbins=nbins, xmin=xmin, xmax=xmax, g_or_l=loc,  level=10)
run_fortran(RUN=RUN_FORT, type="herm", samples=samples, N=N, nbins=nbins, xmin=xmin, xmax=xmax, g_or_l=loc,  level=50)
run_fortran(RUN=RUN_FORT, type="herm", samples=samples, N=N, nbins=nbins, xmin=xmin, xmax=xmax, g_or_l=loc,  level=250)
run_fortran(RUN=RUN_FORT, type="diag", samples=samples, N=N, nbins=nbins, xmin=xmin, xmax=xmax, g_or_l=loc,  level=10)
run_fortran(RUN=RUN_FORT, type="diag", samples=samples, N=N, nbins=nbins, xmin=xmin, xmax=xmax, g_or_l=loc,  level=50)
run_fortran(RUN=RUN_FORT, type="diag", samples=samples, N=N, nbins=nbins, xmin=xmin, xmax=xmax, g_or_l=loc,  level=250)


# plot and fit through gnuplot
run_gnuplot(RUN=RUN_GNU, type="herm", samples=samples, N=N, gnufile="05h.gnu")
run_gnuplot(RUN=RUN_GNU, type="diag", samples=samples, N=N, gnufile="05d.gnu")
