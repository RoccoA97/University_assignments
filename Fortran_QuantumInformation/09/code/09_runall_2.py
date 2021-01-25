import os
import subprocess
import numpy as np



# function for running fortran and converting output into numerical results
def run_fortran(RUN, exec, N, L, K):
    if RUN:
        results = subprocess.run([exec, "-N", str(N),
                                        "-L", str(L),
                                        "-K", str(K)], stdout=subprocess.PIPE)
        results = results.stdout.decode('utf-8').split('\n')[(-K-1):-1]
        results = [float(result.strip()) for result in results]
        return results
    else:
        return


def run_gnuplot(RUN, ofile, K, gnufile):
    cmd  = "gnuplot "
    cmd += "-e \"ofile='" + ofile + "'\" "
    cmd += "-e \"K=" + str(K) + "\" "
    cmd += gnufile
    if RUN:
        os.system(cmd)
    return



RUN_FORT = False
RUN_GNU  = True

comp     = "gfortran"
c_flags  = "-Wall -ffree-line-length-0"
exec     = "./09.o"
src      = "09.f90"
lapack   = "-L/usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0 -llapack"



# compile src
os.system(' '.join([comp, c_flags, "-o", exec, src, lapack]))


# create results directory if it does not exist
if not os.path.isdir('./results'):
    os.system("mkdir results")


# run fortran
N_min  = 3
N_max  = 10
N_step = 1
L_min  = 0.0
L_max  = 3.0
L_step = 0.1
K      = 4

Ns = range(N_min, N_max+1, N_step)
Ls = np.arange(L_min, L_max+L_step, L_step)
Ks = range(K)

Rs = np.zeros((len(Ks),len(Ls),len(Ns)))

if RUN_FORT:
    i = 0
    for L in Ls:
        j = 0
        for N in Ns:
            print("Running for N=", N, "  L=", str(L))
            res = run_fortran(RUN=RUN_FORT, exec=exec, N=N, L=L, K=K)
            for k in Ks:
                Rs[k,i,j] = res[k] / (N-1)
            j += 1
        i += 1

    for k in Ks:
        data = np.column_stack((np.array(Ls),Rs[k]))
        np.savetxt("results/ising_"+str(k)+"_K.dat", data)


# run gnuplot
for k in Ks:
    run_gnuplot(RUN=RUN_GNU, ofile="ising_"+str(k)+"_K.pdf", K=k, gnufile="09_2.gnu")
