import matplotlib.pyplot as plt
import numpy as np
import os
import subprocess



# plot style configuration
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=16)          # controls default text sizes
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=18)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
plt.rc('ytick', labelsize=14)    # fontsize of the tick labels
plt.rc('legend', fontsize=16)    # legend fontsize
plt.rc('figure', titlesize=18)   # fontsize of the figure title



# function for converting fortran output into numerical results
def get_output(cmd, N):
    results = subprocess.run([cmd, str(N)], stdout=subprocess.PIPE)
    results = results.stdout.decode('utf-8').strip().split('\n')
    results = [float(res.strip()) for res in results]
    return np.array(results)





RUN_FORT = False

N_min  = 50
N_max  = 2000
N_step = 30

Ns     = np.int32(np.round(np.logspace(np.log10(N_min), np.log10(N_max), num=N_step, base=10), 0))
Oflags = [" -O0 ", " -O1 ", " -O2 ", " -O3 ", " -Ofast "]
M      = 5

cmd  = "./04.o"
comp = "gfortran"
opt  = " -ffree-line-length-0 "
exec = " -o 04.o "
src  = " 04.f90 "



# get execution times:
if (RUN_FORT):
    matmul_col_data = np.zeros((len(Ns), len(Oflags)*2))
    matmul_row_data = np.zeros((len(Ns), len(Oflags)*2))
    matmul_data     = np.zeros((len(Ns), len(Oflags)*2))

    # No optimization flags
    for j in range(len(Oflags)):
        Oflag = Oflags[j]
        print("["+Oflag+"] : start -------------------------------------------")
        # os.system("rm matmul_col.dat matmul_row.dat matmul.dat")
        os.system(comp + Oflag + opt + exec + src)
        for N,i in zip(Ns,range(len(Ns))):
            print("Running matmul functions for N = " + str(N))
            res = np.zeros((3,M))
            for k in range(M):
                res[:,k] = get_output(cmd, N)
            # store mean
            matmul_col_data[i,j] = np.mean(res[0,:])
            matmul_row_data[i,j] = np.mean(res[1,:])
            matmul_data[i,j]     = np.mean(res[2,:])
            # store standard deviation
            matmul_col_data[i,j+len(Oflags)] = np.std(res[0,:]) / np.sqrt(M)
            matmul_row_data[i,j+len(Oflags)] = np.std(res[1,:]) / np.sqrt(M)
            matmul_data[i,j+len(Oflags)]     = np.std(res[2,:]) / np.sqrt(M)
        print("["+Oflag+"] : end ---------------------------------------------")

    matmul_col_data = np.concatenate((np.array([Ns]).T,matmul_col_data), axis=1)
    matmul_row_data = np.concatenate((np.array([Ns]).T,matmul_row_data), axis=1)
    matmul_data     = np.concatenate((np.array([Ns]).T,matmul_data),     axis=1)

    np.savetxt('matmul_col.dat', matmul_col_data, fmt=' '.join(['%i'] + ['%1.8e']*2*len(Oflags)))
    np.savetxt('matmul_row.dat', matmul_row_data, fmt=' '.join(['%i'] + ['%1.8e']*2*len(Oflags)))
    np.savetxt('matmul.dat',     matmul_data,     fmt=' '.join(['%i'] + ['%1.8e']*2*len(Oflags)))
else:
    matmul_col_data = np.loadtxt("matmul_col.dat")
    matmul_row_data = np.loadtxt("matmul_row.dat")
    matmul_data     = np.loadtxt("matmul.dat")



# Plot the results
for i in range(len(Oflags)):
    Oflag = Oflags[i]
    print("gnuplot -e \"titlename='"+Oflag.strip()+"'\" -e \"mean_col="+str(i+2)+"\" -e \"std_col="+str(i+2+len(Oflags))+"\" 04.p")
    os.system("gnuplot -e \"titlename='"+Oflag.strip()+"'\" -e \"mean_col="+str(i+2)+"\" -e \"std_col="+str(i+2+len(Oflags))+"\" 04.p")
