import matplotlib.pyplot as plt
import numpy as np

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

for order in range(6):
    #Polynomials
    data = np.loadtxt(f"legendre_P{order}_values.txt", skiprows=2)
    ax1.plot(data[:, 0], data[:, 1], label=f'P{order}')
    
    #Derivatives  
    data = np.loadtxt(f"legendre_P{order}_derivatives.txt", skiprows=2)
    ax2.plot(data[:, 0], data[:, 1], label=f"P'{order}")

ax1.set_title('Legendre Polynomials')
ax1.set_ylabel('P_n(x)')
ax1.legend()
ax1.grid(True)

ax2.set_title('Derivatives')
ax2.set_xlabel('x')
ax2.set_ylabel("P'_n(x)")
ax2.legend()
ax2.grid(True)

plt.tight_layout()
plt.savefig('legendre_plot.png', dpi=150)
plt.show()