import matplotlib.pyplot as plt
import numpy as np
import os



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





Ns = [n for n in range(10, 100, 10)] + [n for n in range(100, 1100, 100)]


# get execution times:

# No optimization flags
os.system("rm matmul_col_results.txt matmul_row_results.txt matmul_results.txt")
os.system("gfortran -fno-range-check 01_ex3.f90 -o 01_ex3.o")
for n in Ns:
    os.system("./01_ex3.o " + str(n))

matmul_col_data = np.loadtxt("matmul_col_results.txt")
matmul_row_data = np.loadtxt("matmul_row_results.txt")
matmul_data     = np.loadtxt("matmul_results.txt")


# -O0 optimization flag
os.system("rm matmul_col_results.txt matmul_row_results.txt matmul_results.txt")
os.system("gfortran -fno-range-check -O0 01_ex3.f90 -o 01_ex3.o")
for n in Ns:
    os.system("./01_ex3.o " + str(n))

O0_matmul_col_data = np.loadtxt("matmul_col_results.txt")
O0_matmul_row_data = np.loadtxt("matmul_row_results.txt")
O0_matmul_data     = np.loadtxt("matmul_results.txt")


# -O1 optimization flag
os.system("rm matmul_col_results.txt matmul_row_results.txt matmul_results.txt")
os.system("gfortran -fno-range-check -O1 01_ex3.f90 -o 01_ex3.o")
for n in Ns:
    os.system("./01_ex3.o " + str(n))

O1_matmul_col_data = np.loadtxt("matmul_col_results.txt")
O1_matmul_row_data = np.loadtxt("matmul_row_results.txt")
O1_matmul_data     = np.loadtxt("matmul_results.txt")


# -O2 optimization flag
os.system("rm matmul_col_results.txt matmul_row_results.txt matmul_results.txt")
os.system("gfortran -fno-range-check -O2 01_ex3.f90 -o 01_ex3.o")
for n in Ns:
    os.system("./01_ex3.o " + str(n))

O2_matmul_col_data = np.loadtxt("matmul_col_results.txt")
O2_matmul_row_data = np.loadtxt("matmul_row_results.txt")
O2_matmul_data     = np.loadtxt("matmul_results.txt")


# -O3 optimization flag
os.system("rm matmul_col_results.txt matmul_row_results.txt matmul_results.txt")
os.system("gfortran -fno-range-check -O3 01_ex3.f90 -o 01_ex3.o")
for n in Ns:
    os.system("./01_ex3.o " + str(n))

O3_matmul_col_data = np.loadtxt("matmul_col_results.txt")
O3_matmul_row_data = np.loadtxt("matmul_row_results.txt")
O3_matmul_data     = np.loadtxt("matmul_results.txt")


# -Ofast optimization flag
os.system("rm matmul_col_results.txt matmul_row_results.txt matmul_results.txt")
os.system("gfortran -fno-range-check -Ofast 01_ex3.f90 -o 01_ex3.o")
for n in Ns:
    os.system("./01_ex3.o " + str(n))

Ofast_matmul_col_data = np.loadtxt("matmul_col_results.txt")
Ofast_matmul_row_data = np.loadtxt("matmul_row_results.txt")
Ofast_matmul_data     = np.loadtxt("matmul_results.txt")




# Plot the results

# All
plt.yscale("log")
plt.plot(np.int32(matmul_col_data[:,0]), matmul_col_data[:,1], label="matmul\_col")
plt.plot(np.int32(matmul_row_data[:,0]), matmul_row_data[:,1], label="matmul\_row")
plt.plot(np.int32(matmul_data    [:,0]), matmul_data    [:,1], label="matmul")
plt.xlabel("n")
plt.ylabel("CPU time [s]")
plt.title("Performance comparison")
plt.grid(True, linestyle=':')
plt.legend()
plt.savefig("../report/images/01_ex3_performance.pdf")
plt.show()


# matmul_col
plt.yscale("log")
plt.plot(np.int32(O0_matmul_col_data[:,0]), O0_matmul_col_data[:,1], label="-O0")
plt.plot(np.int32(O1_matmul_col_data[:,0]), O1_matmul_col_data[:,1], label="-O1")
plt.plot(np.int32(O2_matmul_col_data[:,0]), O2_matmul_col_data[:,1], label="-O2")
plt.plot(np.int32(O3_matmul_col_data[:,0]), O3_matmul_col_data[:,1], label="-O3")
plt.plot(np.int32(Ofast_matmul_col_data[:,0]), Ofast_matmul_col_data[:,1], label="-Ofast")
plt.xlabel("n")
plt.ylabel("CPU time [s]")
plt.title("matmul\_col optimization")
plt.grid(True, linestyle=':')
plt.legend()
plt.savefig("../report/images/01_ex3_matmul_col_O.pdf")
plt.show()


# matmul_col
plt.yscale("log")
plt.plot(np.int32(O0_matmul_row_data[:,0]), O0_matmul_row_data[:,1], label="-O0")
plt.plot(np.int32(O1_matmul_row_data[:,0]), O1_matmul_row_data[:,1], label="-O1")
plt.plot(np.int32(O2_matmul_row_data[:,0]), O2_matmul_row_data[:,1], label="-O2")
plt.plot(np.int32(O3_matmul_row_data[:,0]), O3_matmul_row_data[:,1], label="-O3")
plt.plot(np.int32(Ofast_matmul_row_data[:,0]), Ofast_matmul_row_data[:,1], label="-Ofast")
plt.xlabel("n")
plt.ylabel("CPU time [s]")
plt.title("matmul\_row optimization")
plt.grid(True, linestyle=':')
plt.legend()
plt.savefig("../report/images/01_ex3_matmul_row_O.pdf")
plt.show()


# matmul
plt.yscale("log")
plt.plot(np.int32(O0_matmul_data[:,0]), O0_matmul_data[:,1], label="-O0")
plt.plot(np.int32(O1_matmul_data[:,0]), O1_matmul_data[:,1], label="-O1")
plt.plot(np.int32(O2_matmul_data[:,0]), O2_matmul_data[:,1], label="-O2")
plt.plot(np.int32(O3_matmul_data[:,0]), O3_matmul_data[:,1], label="-O3")
plt.plot(np.int32(Ofast_matmul_data[:,0]), Ofast_matmul_data[:,1], label="-Ofast")
plt.xlabel("n")
plt.ylabel("CPU time [s]")
plt.title("matmul optimization")
plt.grid(True, linestyle=':')
plt.legend()
plt.savefig("../report/images/01_ex3_matmul_O.pdf")
plt.show()
