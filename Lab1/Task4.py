import json
import random
import numpy as np
from matplotlib import pyplot as plt

def my_shuffle(arr):
    random.shuffle(arr)
    return arr

def run_simulation(dt, tsteps, kd_prot, kp_mRNA, kd_mRNA, alpha, beta, q_arr, codon_length, Np_0, Nr_0):

    Np_t = np.empty(tsteps, dtype=int)
    Np_t[0] = Np_0

    end_occ = np.zeros(tsteps) 
    
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

        if NmRNA > 0:
            end_occ[i+1] = np.mean(mRNAs[:,codon_length-1])
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

def get_codon_rates(codon_seq, rates_dict):
    codon_length = int(len(codon_seq)/3)

    q_arr = np.empty(codon_length)
    for i in range(codon_length):
        codon = codon_seq[3*i:3*(i+1)]
        q_arr[i] = rates_dict[codon]

    return q_arr

def get_seq_variant(codon_seq, cod_from_amin, amin_from_cod):

    var_seq = ""

    for i in range(0, len(codon_seq), 3):
        codon = codon_seq[i:i+3]
        amin = amin_from_cod[codon]
        possible_codons = cod_from_amin[amin]
        r = random.randint(0, len(possible_codons)-1)
        var_seq += possible_codons[r]
    
    return var_seq



f = open("rates_from_codon.json", "r")
rate_from_codon = json.load(f)
f.close()

f = open("aminoacid_from_codon.json", "r")
amin_from_cod = json.load(f)
f.close()

f = open("codons_from_aminoacid.json", "r")
cod_from_amin = json.load(f)
f.close()


codon_seq = "AGCGCGCGGUCACAACGUUACUGUUAUCGAUCCGGUCGAAAAACUGCUGGCAGUGGGGCAUUACCUCGAAUCUACCGUCGAUAUUGCUGA"

codon_length = int(len(codon_seq)/3)
n_variants = 100

seq_vars = []
q_vars = np.empty((n_variants+1, codon_length))

seq_vars.append(codon_seq)
q_vars[0,:] = get_codon_rates(codon_seq, rate_from_codon)

for i in range(1, n_variants+1):
    seq_vars.append(get_seq_variant(codon_seq, cod_from_amin, amin_from_cod))
    q_vars[i,:] = get_codon_rates(seq_vars[i], rate_from_codon)


kd_prot = 1/1800
kp_mRNA = 1/600
kd_mRNA = 1/300

dt = 0.02
tsteps = 3600 * 3600
Np_0 = 0
Nr_0 = 4*kp_mRNA/kd_mRNA

beta = 2
alpha = 0.9

Np_ts = np.empty((n_variants+1, tsteps))
end_occs = np.empty((n_variants+1, tsteps))

for i, q_var in list(zip(range(n_variants+1), q_vars)):
    Np_ts[i,:], end_occs[i,:] = run_simulation(dt, tsteps, kd_prot, kp_mRNA, kd_mRNA, alpha, beta, q_var, codon_length, Np_0, Nr_0)


np.savetxt("data/q_vars.txt", q_vars)
np.savetxt("data/Np_ts.txt", Np_ts)
np.savetxt("data/end_occs.txt", end_occs)

f_params = open("data/params.txt", "w")
f_params.write(f"alpha = {alpha}\n")
f_params.write(f"beta = {beta}\n")
f_params.write(f"dt = {dt}\n")
f_params.write(f"tsteps = {tsteps}\n")
f_params.write(f"variants = {n_variants}\n")
f_params.close()

