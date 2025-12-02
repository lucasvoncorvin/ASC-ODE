import numpy as np
import matplotlib.pyplot as plt

# Helper to load text files
def load_solution(filename):
    data = np.loadtxt(filename)
    t  = data[:,0]
    x  = data[:,1]
    v  = data[:,2]
    return t, x, v

# Load all RK results
t_RK4,    x_RK4,    v_RK4    = load_solution("output_RK4.txt")
t_Gauss,  x_Gauss,  v_Gauss  = load_solution("output_Gauss4.txt")
t_Radau,  x_Radau,  v_Radau  = load_solution("output_Radau4.txt")

# TIME EVOLUTION OF POSITION x(t)
plt.figure(figsize=(10,6))
plt.plot(t_RK4, x_RK4, 'o-', label="RK4 x(t)")
plt.plot(t_Gauss, x_Gauss, 's-', label="Gauss(4) x(t)")
plt.plot(t_Radau, x_Radau, 'd-', label="Radau(4) x(t)")
plt.xlabel("t")
plt.ylabel("x(t)")
plt.title("Time Evolution of Position for N = 10")
plt.grid(True)
plt.legend()

plt.savefig("time_evolution_position_10.png", dpi=300)
print("Saved time_evolution_position_10.png")

# ---------------------------------------------------------
# TIME EVOLUTION OF VELOCITY v(t)
# ---------------------------------------------------------
plt.figure(figsize=(10,6))
plt.plot(t_RK4, v_RK4, 'o-', label="RK4 v(t)")
plt.plot(t_Gauss, v_Gauss, 's-', label="Gauss(4) v(t)")
plt.plot(t_Radau, v_Radau, 'd-', label="Radau(4) v(t)")
plt.xlabel("t")
plt.ylabel("v(t)")
plt.title("Time Evolution of Velocity")
plt.grid(True)
plt.legend()

plt.savefig("time_evolution_velocity_10.png", dpi=300)
print("Saved time_evolution_velocity_10.png")

# ---------------------------------------------------------
# PHASE PLOT x(t) vs v(t)
# ---------------------------------------------------------
plt.figure(figsize=(8,8))
plt.plot(x_RK4,    v_RK4,    'o-', label="RK4")
plt.plot(x_Gauss,  v_Gauss,  's-', label="Gauss(4)")
plt.plot(x_Radau,  v_Radau,  'd-', label="Radau(4)")
plt.xlabel("x(t)")
plt.ylabel("v(t)")
plt.title("Phase Plot N = 10")
plt.grid(True)
plt.legend()

plt.savefig("phase_plot_10.png", dpi=300)
print("Saved phase_plot_10.png")

plt.show()
