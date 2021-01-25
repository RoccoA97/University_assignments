import os
import subprocess
import numpy as np



# function for running fortran and converting output into numerical results
def run_fortran(RUN, exec, N, D, S, M):
    if RUN:
        results = subprocess.run([exec, "-N", str(N),
                                        "-D", str(D),
                                        "-S", str(S),
                                        "-M", str(M)], stdout=subprocess.PIPE)
        results = results.stdout.decode('utf-8').split('\n')[-2].strip()
        results = float(results)
        return results
    else:
        return


def run_gnuplot(RUN, ofile, gnufile):
    cmd  = "gnuplot -e "
    cmd += "\"ofile='" + ofile + "'\" "
    cmd += gnufile
    if RUN:
        os.system(cmd)
    return



RUN_FORT = True
RUN_GNU  = True

comp     = "gfortran"
c_flags  = "-Wall -ffree-line-length-0"
exec     = "./08.o"
src      = "08.f90"



# compile src
os.system(' '.join([comp, c_flags, "-o", exec, src]))


# create results directory if it does not exist
if not os.path.isdir('./results'):
    os.system("mkdir results")


# separable case
N_min  = 10000
N_max  = 1000000
N_step = 10000
D_min  = 2
D_max  = 8
D_step = 2
sep    = 1
mode   = 1

for D in range(D_min, D_max+1, D_step):
    Ns = []
    Ts = []
    for N in range(N_min, N_max+1, N_step):
        print("Running for separable state with N=", N, "  D=", str(D))
        Ns.append(N)
        Ts.append(run_fortran(RUN=RUN_FORT, exec=exec, N=N, D=D, S=sep, M=mode))
    np.savetxt("results/separable_"+str(D)+"_D.dat", np.array([Ns,Ts]).T, fmt=' '.join(['%i'] + ['%1.8e']))


# non-separable case
N_min  = 2
N_max  = 12
N_step = 1
D_min  = 2
D_max  = 5
D_step = 1
sep    = 0
mode   = 1

for D in range(D_min, D_max+1, D_step):
    Ns = []
    Ts = []
    for N in range(N_min, N_max+1, N_step):
        print("Running for non-separable state with N=", N, "  D=", str(D))
        Ns.append(N)
        Ts.append(run_fortran(RUN=RUN_FORT, exec=exec, N=N, D=D, S=sep, M=mode))
    np.savetxt("results/nonseparable_"+str(D)+"_D.dat", np.array([Ns,Ts]).T, fmt=' '.join(['%i'] + ['%1.8e']))


# non-separable case, d=2 (more samples)
N_min  = 2
N_max  = 28
N_step = 1
D_min  = 2
D_max  = 2
D_step = 1
sep    = 0
mode   = 1

for D in range(D_min, D_max+1, D_step):
    Ns = []
    Ts = []
    for N in range(N_min, N_max+1, N_step):
        print("Running for non-separable state with N=", N, "  D=", str(D))
        Ns.append(N)
        Ts.append(run_fortran(RUN=RUN_FORT, exec=exec, N=N, D=D, S=sep, M=mode))
    np.savetxt("results/nonseparable_"+str(D)+"_D_2.dat", np.array([Ns,Ts]).T, fmt=' '.join(['%i'] + ['%1.8e']))



# run gnuplot
run_gnuplot(RUN=RUN_GNU, ofile="separable.pdf",      gnufile="08_s.gnu"   )
run_gnuplot(RUN=RUN_GNU, ofile="nonseparable.pdf",   gnufile="08_ns.gnu"  )
run_gnuplot(RUN=RUN_GNU, ofile="nonseparable_2.pdf", gnufile="08_ns_2.gnu")
