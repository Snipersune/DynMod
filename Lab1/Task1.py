import random
import numpy as np
from matplotlib import pyplot as plt


def run_simulation(dt, tsteps, kp, kd, proteins_per_mRNA, Np_0):
    Np_t = np.empty(tsteps, dtype=int)
    Np_t[0] = Np_0 

    Np = Np_0
    for i in range(tsteps-1):
        r = random.random()
        if r < dt*kd*Np:
            Np -= 1

        r = random.random()
        if r < dt*kp:
            Np += proteins_per_mRNA
        
        Np_t[i+1] = Np

    return Np_t

kp = 1/600
kd = 1/1800
Np_0 = 0

dt = 1
tsteps = 500 * 3600


prot_per_mRNA = 1
Np_t1 = run_simulation(dt, tsteps, kp, kd, prot_per_mRNA, Np_0)

prot_per_mRNA = 10
Np_t10 = run_simulation(dt, tsteps, kp, kd, prot_per_mRNA, Np_0)

Np1_mean = np.mean(Np_t1)
Np1_var = np.var(Np_t1)
ff1 = Np1_var/Np1_mean

Np10_mean = np.mean(Np_t10)
Np10_var = np.var(Np_t10)
ff10 = Np10_var/Np10_mean

print("Proteins per mRNA = 1")
print("----------------------")
print("Mean:", Np1_mean)
print("Variance:", Np1_var)
print("Fano factor:", ff1)
print("")
print("Proteins per mRNA = 10")
print("----------------------")
print("Mean:", Np10_mean)
print("Variance:", Np10_var)
print("Fano factor:", ff10)
print("")

print(ff10/prot_per_mRNA)

plt.figure()
plt.plot(range(0,tsteps*dt, dt), Np_t1, label="prot_per_mRNA=1")
plt.plot(range(0,tsteps*dt, dt), Np_t10, label="prot_per_mRNA=10")
plt.legend()
plt.show()

