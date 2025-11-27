import numpy as np
import matplotlib.pyplot as plt

# LOAD DATA
exp  = np.loadtxt("rc_explicit.txt")
imp  = np.loadtxt("rc_implicit.txt")
cran = np.loadtxt("rc_crank.txt")

# Columns: time, UC
t_exp,  U_exp  = exp[:,0], exp[:,1]
t_imp,  U_imp  = imp[:,0], imp[:,1]
t_cran, U_cran = cran[:,0], cran[:,1]


# COMPUTE NUMERICAL DERIVATIVES FOR PHASE PLOTS
def numerical_derivative(t, y):
    dy = np.zeros_like(y)
    dy[1:-1] = (y[2:] - y[:-2]) / (t[2:] - t[:-2])     # centered difference
    dy[0]    = (y[1] - y[0]) / (t[1] - t[0])
    dy[-1]   = (y[-1] - y[-2]) / (t[-1] - t[-2])
    return dy

dU_exp  = numerical_derivative(t_exp,  U_exp)
dU_imp  = numerical_derivative(t_imp,  U_imp)
dU_cran = numerical_derivative(t_cran, U_cran)


# FIGURE 1: TIME PLOT
plt.figure(figsize=(10,6))
#plt.plot(t_exp,  U_exp,  label="Explicit Euler")
plt.plot(t_imp,  U_imp,  label="Implicit Euler")
plt.plot(t_cran, U_cran, label="Crank–Nicolson")

plt.xlabel("time t")
plt.ylabel("Voltage UC(t)")
plt.title("RC Circuit – UC(t) for Three Methods")

plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("rc_methods_time_1_50_R_C_stiff_woEE.png", dpi=200)
plt.show()


# FIGURE 2: PHASE PLOT
plt.figure(figsize=(10,6))

#plt.plot(U_exp,  dU_exp,  label="Explicit Euler")
plt.plot(U_imp,  dU_imp,  label="Implicit Euler")
plt.plot(U_cran, dU_cran, label="Crank–Nicolson")

plt.xlabel("UC(t)")
plt.ylabel("dUC/dt")
plt.title("Phase Plot of RC Circuit")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("rc_methods_phase_1_50_R_C_stiff.png", dpi=200)
plt.show()
