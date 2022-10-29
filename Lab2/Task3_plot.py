from cProfile import label
import numpy as np
from matplotlib import pyplot as plt


gamma = np.loadtxt("data/task3_gammas.txt")
nMU = []
for i in range(len(gamma)):
    filename = f"data/task3_nMUg{i+1}.txt"
    nMU.append(np.loadtxt(filename))
    
print(len(nMU))

plt.figure()
plt.plot(nMU[0][0,:], label="U")
plt.plot(nMU[0][1,:], label="M")
plt.legend()

plt.figure()
plt.plot(nMU[4][0,:], label="U")
plt.plot(nMU[4][1,:], label="M")
plt.legend()

diffMA1 = 2*nMU[0][1,:] + nMU[0][0,:] - 60
diffMA2 = 2*nMU[0][1,:] + nMU[0][0,:] - 60

for i in range(len(nMU)):
    diffMA = 2*nMU[i][1,:] + nMU[i][0,:] - 60

    plt.figure()
    plt.hist(diffMA, bins=120, label=f"$\gamma = {gamma[i]}$")
    plt.legend()

plt.show()
