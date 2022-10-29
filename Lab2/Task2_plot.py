from cProfile import label
import numpy as np
from matplotlib import pyplot as plt


nMU1 = np.loadtxt("data/task2.1_F10.txt")
nMU2 = np.loadtxt("data/task2.2_F10.txt")


plt.figure()
plt.plot(nMU1[0,:], label="U")
plt.plot(nMU1[1,:], label="M")
plt.legend()

plt.figure()
plt.plot(nMU2[0,:], label="U")
plt.plot(nMU2[1,:], label="M")
plt.legend()

diffMA1 = 2*nMU1[1,:] + nMU1[0,:] - 60
diffMA2 = 2*nMU2[1,:] + nMU2[0,:] - 60

plt.figure()
plt.hist(diffMA1, bins=120, label="F=2")
plt.legend()

plt.figure()
plt.hist(diffMA2, bins=120, label="F=4")
plt.legend()

plt.show()