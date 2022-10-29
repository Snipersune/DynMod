import code
import random
import numpy as np
from matplotlib import pyplot as plt


def run_simulation(dt, tsteps, kd, alpha, beta, q_arr, codon_length, Np_0):
    mRNA_occupancy = np.zeros(codon_length, dtype=int)
    end_occ = np.empty(tsteps, dtype=int) 
    end_occ[0] = 0

    Np_t = np.empty(tsteps, dtype=int)
    Np_t[0] = Np_0

    Np = Np_0
    for i in range(tsteps-1):

        r_decay = random.random()
        if r_decay < dt*kd*Np:
            Np -= 1
        
        mRNA_occupancy, Np = update_occupancy(mRNA_occupancy, Np, dt, alpha, beta, q_arr, codon_length)

        end_occ[i+1] = mRNA_occupancy[codon_length-1]
        Np_t[i+1] = Np

    return Np_t, end_occ

def update_occupancy(mRNA, Np, dt, alpha, beta, q_arr, codon_length):

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
        if r_attach < alpha*dt:
            mRNA[0] = 1
    
    if end_occupied == 1:
        r_detach = random.random()
        if r_detach < beta*dt:
            mRNA[codon_length-1] = 0
            Np += 1
    
    return mRNA, Np
    

kd = 1/1800
codon_length = 60
q_arr = np.ones(codon_length)

dt = 0.25
tsteps = 48 * 3600
Np_0 = 0

beta = [0.5, 0.25]
alpha = np.arange(0.1, 1.0, 0.1)

Np_t05 = np.empty((tsteps, len(alpha)), dtype=int)
end_occ05 = np.empty((tsteps, len(alpha)), dtype=int)

Np_t025 = np.empty((tsteps, len(alpha)), dtype=int)
end_occ025 = np.empty((tsteps, len(alpha)), dtype=int)

for i in range(len(alpha)):
    Np_t05[:,i], end_occ05[:,i] = run_simulation(dt, tsteps, kd, alpha[i], beta[0], q_arr, codon_length, Np_0)

for i in range(len(alpha)):
    Np_t025[:,i], end_occ025[:,i] = run_simulation(dt, tsteps, kd, alpha[i], beta[1], q_arr, codon_length, Np_0)


tsteps_transient = int(tsteps/5)

Np_t05_noTransient = Np_t05[tsteps_transient:,:]
Np_t025_noTransient = Np_t025[tsteps_transient:,:]

end_occ05 = end_occ05[tsteps_transient:,:]
end_occ025 = end_occ025[tsteps_transient:,:]

pL05 = np.mean(end_occ05, axis = 0)
pL025 = np.mean(end_occ025, axis = 0)

current05 = beta[0]*pL05
current025 = beta[1]*pL025

mean05 = np.mean(Np_t05_noTransient, axis=0)
variance05 = np.var(Np_t05_noTransient, axis=0)

mean025 = np.mean(Np_t025_noTransient, axis=0)
variance025 = np.var(Np_t025_noTransient, axis=0)

ff05 = variance05/mean05
ff025 = variance025/mean025

plt.figure()
plt.plot(np.arange(0,tsteps*dt, dt), Np_t05)

plt.figure()
plt.plot(alpha, current05)
plt.plot(alpha, current025)

plt.show()