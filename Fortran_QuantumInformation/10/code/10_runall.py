import os
import subprocess
import numpy as np



# function for running fortran and converting output into numerical results
def run_fortran(RUN, exec, N, L, I):
    if RUN:
        results = subprocess.run([exec, "-N", str(N),
                                        "-L", str(L),
                                        "-I", str(I)], stdout=subprocess.PIPE)
        results = results.stdout.decode('utf-8').split('\n')[-2]
        results = float(results.strip())
        return results
    else:
        return


def run_gnuplot(RUN, ofile, I, N_min, N_max, gnufile):
    cmd  = "gnuplot "
    cmd += "-e \"ofile='" + ofile + "'\" "
    cmd += "-e \"N_min=" + str(N_min) + "\" "
    cmd += "-e \"N_max=" + str(N_max) + "\" "
    cmd += "-e \"I=" + str(I) + "\" "
    cmd += gnufile
    if RUN:
        os.system(cmd)
    return



RUN_FORT = False
RUN_GNU  = True

comp     = "gfortran"
c_flags  = "-Wall -ffree-line-length-0"
exec     = "./10.o"
src      = "10.f90"
lapack   = "-L/usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0 -llapack"



# compile src
os.system(' '.join([comp, c_flags, "-o", exec, src, lapack]))


# create results directory if it does not exist
if not os.path.isdir('./results'):
    os.system("mkdir results")


# run fortran
N_min  = 2
N_max  = 4
N_step = 1
L_min  = -3.0
L_max  = 3.0
L_step = 0.1
I      = 100

Ns = range(N_min, N_max+1, N_step)
Ls = np.round(np.arange(L_min, L_max+L_step, L_step), 2)

Rs = np.zeros((len(Ls),len(Ns)))

if RUN_FORT:
    j = 0
    for N in Ns:
        i = 0
        for L in Ls:
            print("Running for N=", N, "  L=", str(L))
            res = run_fortran(RUN=RUN_FORT, exec=exec, N=N, L=L, I=I)
            Rs[i,j] = res
            i += 1
        j += 1

    data = np.column_stack((np.array(Ls),Rs))
    np.savetxt("results/ising_RSRG_"+str(I)+"_I.dat", data)


# run gnuplot
run_gnuplot(RUN=RUN_GNU, ofile="ising_RSRG_"+str(I)+"_I.pdf", N_min=N_min, N_max=N_max, I=I, gnufile="10.gnu")
