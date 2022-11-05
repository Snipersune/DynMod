import numpy as np
import matplotlib.pyplot as plt


# Load data from files
rand_data = np.loadtxt("data/rand_frac.txt")
sorted_data = np.loadtxt("data/sorted_frac.txt")

plt.figure()
plt.plot(rand_data[0,:], rand_data[1,:], label="ecoli")
plt.plot(rand_data[0,:], rand_data[2,:], label="kidney")
plt.plot(rand_data[0,:], rand_data[3,:], label="yeast")
plt.xlabel("f")
plt.ylabel("D")
plt.legend()

plt.figure()
plt.plot(sorted_data[0,:], sorted_data[1,:], label="ecoli")
plt.plot(sorted_data[0,:], sorted_data[2,:], label="kidney")
plt.plot(sorted_data[0,:], sorted_data[3,:], label="yeast")
plt.xlabel("f")
plt.ylabel("D")
plt.legend()
plt.show()

plt.figure()
plt.plot(rand_data[0,:], rand_data[1,:], label="random")
plt.plot(sorted_data[0,:], sorted_data[1,:], label="highest degrees")
plt.xlabel("f")
plt.ylabel("D")
plt.legend()

plt.figure()
plt.plot(rand_data[0,:], rand_data[2,:], label="random")
plt.plot(sorted_data[0,:], sorted_data[2,:], label="highest degrees")
plt.xlabel("f")
plt.ylabel("D")
plt.legend()

plt.figure()
plt.plot(rand_data[0,:], rand_data[3,:], label="random")
plt.plot(sorted_data[0,:], sorted_data[3,:], label="highest degrees")
plt.xlabel("f")
plt.ylabel("D")
plt.legend()
plt.show()


