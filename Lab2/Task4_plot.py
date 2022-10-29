from cProfile import label
import numpy as np
from matplotlib import pyplot as plt


div_TS = np.loadtxt("data/task4_divTS.txt")
nMU = np.loadtxt("data/task4_F6.txt")

plt.figure()
plt.plot(nMU[0,:], label="U")
plt.plot(nMU[1,:], label="M")
plt.vlines(div_TS, ymin=0, ymax=40, colors="red", linestyles="dashed")
plt.legend()

diffMA = 2*nMU[1,:] + nMU[0,:] - 60

plt.figure()
plt.hist(diffMA, bins=120, label=f"$\gamma = -1.0$")
plt.legend()

plt.show()
