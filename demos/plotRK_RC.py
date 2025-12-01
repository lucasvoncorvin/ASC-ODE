import numpy as np
import matplotlib.pyplot as plt

 
# Load the solution file
def load_solution(filename):
    data = np.loadtxt(filename)
    t  = data[:,0]
    UC = data[:,1]
    return t, UC

# Numerical derivative (central finite differences)
def numerical_derivative(t, y):
    dy = np.zeros_like(y)
    dt = np.diff(t)

    # central differences
    dy[1:-1] = (y[2:] - y[:-2]) / (t[2:] - t[:-2])

    # forward/backward for edges
    dy[0]  = (y[1] - y[0]) / (t[1] - t[0])
    dy[-1] = (y[-1] - y[-2]) / (t[-1] - t[-2])

    return dy

# Load data
t_RK4,   UC_RK4   = load_solution("output_RK4.txt")
t_Gauss, UC_Gauss = load_solution("output_Gauss4.txt")
t_Radau, UC_Radau = load_solution("output_Radau4.txt")

# Numerical derivatives
dU_RK4   = numerical_derivative(t_RK4,   UC_RK4)
dU_Gauss = numerical_derivative(t_Gauss, UC_Gauss)
dU_Radau = numerical_derivative(t_Radau, UC_Radau)

plt.figure(figsize=(10,6))

#plt.plot(t_RK4,   UC_RK4,   label="RK4",      marker="o", markersize=3)
plt.plot(t_Gauss, UC_Gauss, label="Gauss(4)", marker="s", markersize=3)
plt.plot(t_Radau, UC_Radau, label="Radau(4)", marker="d", markersize=3)

plt.xlabel("Time t")
plt.ylabel("Capacitor Voltage U_C(t)")
plt.title("RC Circuit: Time Evolution of U_C N =1000")
plt.grid(True)
plt.legend()

plt.savefig("UC_time_evolution_stiff_1000_NORK4.png", dpi=300)
print("Saved UC_time_evolution.png")

# Phase plot: U_C vs dU_C/dt

plt.figure(figsize=(8,6))
#plt.plot(UC_RK4,   dU_RK4,   label="RK4", marker="o", markersize=3)
plt.plot(UC_Gauss, dU_Gauss, label="Gauss(4)", marker="s", markersize=3)
plt.plot(UC_Radau, dU_Radau, label="Radau(4)", marker="d", markersize=3)

plt.xlabel("U_C (capacitor voltage)")
plt.ylabel("dU_C/dt (numerical derivative)")
plt.title("Phase Plot: U_C vs dU_C/dt N = 1000")
plt.grid(True)
plt.legend()

plt.savefig("phase_plot_U_dU_stiff_1000_NORK4.png", dpi=300)
print("Saved phase_plot_U_dU.png")

plt.show()
