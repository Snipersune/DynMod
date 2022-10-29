import code
import random
import numpy as np
from matplotlib import pyplot as plt

def my_shuffle(arr):
    random.shuffle(arr)
    return arr

def run_simulation(dt, tsteps, kd_prot, kp_mRNA, kd_mRNA, alpha, beta, q_arr, codon_length, Np_0, Nr_0):

    Np_t = np.empty(tsteps, dtype=int)
    Np_t[0] = Np_0

    sum = 0
    Np = Np_0
    Nr = Nr_0
    NmRNA = 0
    mRNAs = np.empty(0)
    for i in range(tsteps-1):

        # Simulate mRNA decay
        r_kd_mRNA = random.random()
        if r_kd_mRNA < kd_mRNA*NmRNA*dt:
            remove_idx = random.randint(0,NmRNA-1)  # mRNA to remove
            Nr += np.sum(mRNAs[remove_idx,:])       # Add back ribosomes to free bulk
            mRNAs = np.delete(mRNAs, remove_idx, axis=0)
            NmRNA -= 1

        # Simulate mRNA production
        r_kp_mRNA = random.random()
        if r_kp_mRNA < kp_mRNA*dt:
            if NmRNA == 0:
                mRNAs = np.append(mRNAs, np.zeros(codon_length))
                mRNAs = np.reshape(mRNAs, (1,-1))
            else:
                mRNAs = np.vstack((mRNAs, np.zeros(codon_length)))
            NmRNA += 1
        
        # Simulate protein decay
        r_decay = random.random()
        if r_decay < dt*kd_prot*Np:
            Np -= 1
        
        # Simulate ribosome dynamics
        for idx in my_shuffle(np.arange(NmRNA)):    # Randomize order in which to check mRNAs 
            mRNAs[idx,:], Np, Nr = update_occupancy(mRNAs[idx,:], Np, Nr, dt, alpha, beta, q_arr, codon_length)

        Np_t[i+1] = Np

    return Np_t

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
    

kd_prot = 1/1800
kp_mRNA = 1/600
kd_mRNA = 1/300
codon_length = 60
q_arr = np.ones(codon_length)

dt = 0.1
tsteps = 720 * 3600
Np_0 = 0
Nr_0 = 4*kp_mRNA/kd_mRNA

beta = 2
alpha = 0.9

Np_t = run_simulation(dt, tsteps, kd_prot, kp_mRNA, kd_mRNA, alpha, beta, q_arr, codon_length, Np_0, Nr_0)

tsteps_transient = int(tsteps/10)

Np_t_noTransient = Np_t[tsteps_transient:]

mean = np.mean(Np_t_noTransient)
variance = np.var(Np_t_noTransient)
ff = variance/mean

print("Fano factor:", ff)

plt.figure()
plt.plot(np.arange(0,tsteps*dt, dt), Np_t)
plt.show()