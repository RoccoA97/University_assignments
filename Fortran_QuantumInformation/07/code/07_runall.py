import os
import subprocess



def run_fortran(RUN, NX, L, NT, T, TF, M, O):
    if RUN:
        subprocess.run([exec, "-NX", str(NX),
                              "-L",  str(L),
                              "-NT", str(NT),
                              "-T",  str(T),
                              "-TF", str(TF),
                              "-M",  str(M),
                              "-O",  str(O)])
    return

def run_gnuplot(RUN, ofile, gnufile):
    print("\nNow fitting and plotting:")
    cmd  = "gnuplot -e "
    cmd += "\"ofile='" + ofile + "'\" "
    cmd += gnufile
    if RUN:
        os.system(cmd)
    return



comp     = "gfortran"
c_flags  = "-Wall -ffree-line-length-0"
exec     = "./07.o"
src      = "07.f90"
lapack   = "-L/usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0 -llapack"
fftw3    = "-L/usr/lib/x86_64-linux-gnu -lfftw3 -lfftw3f"

NX = 2000
L  =   20
NT = 2000
T  =    1
TF =    8



# compile src
os.system(' '.join([comp, c_flags, "-o", exec, src, lapack, fftw3]))



# run executable
M = 1
O = 1
run_fortran(RUN=True, NX=NX, L=L, NT=NT, T=T, TF=TF, M=M, O=O)
run_gnuplot(RUN=True, ofile="cmap_"+str(M)+"_M_"+str(O)+"_O.pdf", gnufile="07_cmap.gnu")
run_gnuplot(RUN=True, ofile="expv_"+str(M)+"_M_"+str(O)+"_O.pdf", gnufile="07_expv.gnu")
run_gnuplot(RUN=True, ofile="gif_" +str(M)+"_M_"+str(O)+"_O.gif", gnufile="07_gif.gnu" )

M = 1
O = 5
run_fortran(RUN=True, NX=NX, L=L, NT=NT, T=T, TF=TF, M=M, O=O)
run_gnuplot(RUN=True, ofile="cmap_"+str(M)+"_M_"+str(O)+"_O.pdf", gnufile="07_cmap.gnu")
run_gnuplot(RUN=True, ofile="expv_"+str(M)+"_M_"+str(O)+"_O.pdf", gnufile="07_expv.gnu")
run_gnuplot(RUN=True, ofile="gif_" +str(M)+"_M_"+str(O)+"_O.gif", gnufile="07_gif.gnu" )

M = 5
O = 1
run_fortran(RUN=True, NX=NX, L=L, NT=NT, T=T, TF=TF, M=M, O=O)
run_gnuplot(RUN=True, ofile="cmap_"+str(M)+"_M_"+str(O)+"_O.pdf", gnufile="07_cmap.gnu")
run_gnuplot(RUN=True, ofile="expv_"+str(M)+"_M_"+str(O)+"_O.pdf", gnufile="07_expv.gnu")
run_gnuplot(RUN=True, ofile="gif_" +str(M)+"_M_"+str(O)+"_O.gif", gnufile="07_gif.gnu" )

M = 0.2
O = 5
run_fortran(RUN=True, NX=NX, L=L, NT=NT, T=T, TF=TF, M=M, O=O)
run_gnuplot(RUN=True, ofile="cmap_"+str(M)+"_M_"+str(O)+"_O.pdf", gnufile="07_cmap.gnu")
run_gnuplot(RUN=True, ofile="expv_"+str(M)+"_M_"+str(O)+"_O.pdf", gnufile="07_expv.gnu")
run_gnuplot(RUN=True, ofile="gif_" +str(M)+"_M_"+str(O)+"_O.gif", gnufile="07_gif.gnu" )
