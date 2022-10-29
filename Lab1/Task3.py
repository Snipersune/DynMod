import code
import random
import numpy as np
from matplotlib import pyplot as plt


def run_simulation(dt, tsteps, kd, alpha, beta, q_arr, codon_length, Np_0, Nr_0):
    mRNA_occupancy = np.zeros(codon_length, dtype=int)
    end_occ = np.empty(tsteps, dtype=int) 
    end_occ[0] = 0

    Np_t = np.empty(tsteps, dtype=int)
    Np_t[0] = Np_0

    Np = Np_0
    Nr = Nr_0
    for i in range(tsteps-1):

        r_decay = random.random()
        if r_decay < dt*kd*Np:
            Np -= 1
        
        mRNA_occupancy, Np, Nr = update_occupancy(mRNA_occupancy, Np, Nr, dt, alpha, beta, q_arr, codon_length)

        end_occ[i+1] = mRNA_occupancy[codon_length-1]
        Np_t[i+1] = Np

    return Np_t, end_occ

def update_occupancy(mRNA, Np, Nr, dt, alpha, beta, q_arr, codon_length):

    # Check initial occupancy at start and end before any movement
    start_occupied = mRNA[0]
    end_occupied = mRNA[codon_length-1]

    # Find indeces of occupied codons except end site and 
    # shuffle list for simulation purposes
    occupied_idx = np.where(mRNA[0:codon_length-1] == 1)[0]
    random.shuffle(occupied_idx)

    # Simulate movement
    for i in occupied_idx:
        if mRNA[i+1] == 0:
            r_step = random.random()
            if r_step < q_arr[i]*dt:
                mRNA[i] = 0
                mRNA[i+1] = 1

    # Allow for attachment and detachment depending on initial occupancy
    if start_occupied == 0:
        r_attach = random.random()
        if r_attach < alpha*Nr*dt:
            mRNA[0] = 1
            Nr -= 1
    
    if end_occupied == 1:
        r_detach = random.random()
        if r_detach < beta*dt:
            mRNA[codon_length-1] = 0
            Np += 1
            Nr += 1
    
    return mRNA, Np, Nr
    

kd = 1/1800
codon_length = 60
q_arr = np.ones(codon_length)

dt = 0.1
tsteps = 720 * 3600
Np_0 = 0
Nr_0 = 4

beta = 2
alpha = 0.9

Np_t, end_occ = run_simulation(dt, tsteps, kd, alpha, beta, q_arr, codon_length, Np_0, Nr_0)

tsteps_transient = int(tsteps/10)

Np_t_noTransient = Np_t[tsteps_transient:]
end_occ = end_occ[tsteps_transient:]

pL05 = np.mean(end_occ)
current05 = beta*pL05

mean = np.mean(Np_t_noTransient)
variance = np.var(Np_t_noTransient)
ff = variance/mean

print("Fano factor:", ff)

plt.figure()
plt.plot(np.arange(0,tsteps*dt, dt), Np_t)
plt.show()