from cProfile import label
import numpy as np
from matplotlib import pyplot as plt


nMUF2 = np.loadtxt("data/task1_F2.txt")
nMUF4 = np.loadtxt("data/task1_F4.txt")
nMUF6 = np.loadtxt("data/task1_F6.txt")

plt.figure()
plt.plot(nMUF2[0,:], label="U")
plt.plot(nMUF2[1,:], label="M")
plt.legend()

plt.figure()
plt.plot(nMUF4[0,:], label="U")
plt.plot(nMUF4[1,:], label="M")
plt.legend()

plt.figure()
plt.plot(nMUF6[0,:], label="U")
plt.plot(nMUF6[1,:], label="M")
plt.legend()

diffMA2 = 2*nMUF2[1,:] + nMUF2[0,:] - 60
diffMA4 = 2*nMUF4[1,:] + nMUF4[0,:] - 60
diffMA6 = 2*nMUF6[1,:] + nMUF6[0,:] - 60 

plt.figure()
plt.hist(diffMA2, bins=120, label="F=2")
plt.legend()

plt.figure()
plt.hist(diffMA4, bins=120, label="F=4")
plt.legend()

plt.figure()
plt.hist(diffMA6, bins=120, label="F=6")
plt.legend()

plt.show()