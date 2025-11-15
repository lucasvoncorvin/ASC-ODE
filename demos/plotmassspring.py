import numpy as np
data = np.loadtxt('output_test_ode.txt', usecols=(0, 1, 2))
# print (data)

import matplotlib.pyplot as plt

# Figure 1
plt.figure(figsize=(7,5))
plt.plot(data[:,0], data[:,1], label='position')
plt.plot(data[:,0], data[:,2], label='velocity')
plt.xlabel('time')
plt.ylabel('value')
plt.title('Mass-Spring System Time Evolution')
plt.legend()
plt.grid()

plt.savefig("time_evolution_100_CN.png", dpi=300, bbox_inches='tight')


# Figure 2 
plt.figure(figsize=(7,5))
plt.plot(data[:,1], data[:,2], label='phase plot')
plt.xlabel('position')
plt.ylabel('velocity')
plt.title('Mass-Spring System Phase Plot')
plt.legend()
plt.grid()

plt.savefig("phase_plot_100_CN.png", dpi=300, bbox_inches='tight')

plt.show()