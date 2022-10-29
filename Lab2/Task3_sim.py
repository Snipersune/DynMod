import numpy as np
import random
import math 

def fromFToAlpha(F):
    return F/(1.+F)

def getRandState(L):
    state = np.empty(L, dtype=int)
    for i in range(L):
        state[i] = random.randint(-1,1)
    return state

def getPc(L, gamma):
    Pc = np.empty(L-1)
    Pc[0] = 1
    for i in range(1, L-1):
        Pc[i] = Pc[i-1] + (i+1)**gamma
    return Pc

def localInteraction(Pc, n1, L):

    sym_dist = min(n1, L-1-n1)
    max_dist = max(n1, L-1-n1)

    P_tot = Pc[max_dist-1] + Pc[sym_dist-1]
    
    r_dist = P_tot*random.random()

    n2 = -1
    d = 1
    while d <= sym_dist:
        if r_dist < 2*Pc[d-1]:
            n2 = n1 + d*(1 - 2*random.randint(0,1))
            break
        d += 1
    
    if n2 == -1:
        while d <= max_dist:
            if r_dist < Pc[sym_dist-1] + Pc[d-1]:
                if n1 <= max_dist:
                    n2 = n1 - d
                else:
                    n2 = n1 + d
                break
            d += 1
    
    return n2

def performStep(state, L, alpha, Pc):

    for i in range(L):

        # Feedback effect
        r_feedback = random.random()
        if r_feedback < alpha:
            n1 = random.randint(0, L-1)
            n2 = localInteraction(Pc, n1, L)

            if state[n1] != state[n2]:
                state[n1] += state[n2]
        
        # Noice effect
        else:
            n1 = random.randint(0, L-1)
            if state[n1] != 0:
                state[n1] = 0
            else:
                r_A = random.random()
                if r_A < 0.5:
                    state[n1] = -1
                else:
                    state[n1] = 1

    return state
        
random.seed(22)

L = 60
F = 6
gamma = np.linspace(-3.0, -1.0, num=5)

tsteps = 100000

# To store number of M/A and U at each timestep
nMU = np.empty((len(gamma), 2, tsteps), dtype=int)

alpha = fromFToAlpha(F)

for idx, gma in list(zip(range(len(gamma)), gamma)):
    state = getRandState(L)
    Pc = getPc(L, gma)
    for t in range(tsteps):
        state = performStep(state, L, alpha, Pc)
        nMU[idx,0,t] = sum(state == 0)    
        nMU[idx,1,t] = sum(state == 1)

np.savetxt("data/task3_gammas.txt", gamma)
for i in range(len(gamma)):
    filename = f"data/task3_nMUg{i+1}.txt"
    np.savetxt(filename, nMU[i,:,:], fmt="%d")
