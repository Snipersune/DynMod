import numpy as np
import random


def fromFToAlpha(F):
    return F/(1.+F)

def getRandState(L):
    state = np.empty(L, dtype=int)
    for i in range(L):
        state[i] = random.randint(-1,1)
    return state

def performStep1(state, L, alpha):

    for i in range(L):

        # Feedback effect
        r_feedback = random.random()
        if r_feedback < alpha:
            n1 = random.randint(0, L-1)
            if state[n1] == 0:
                n2 = random.randint(0, L-1)
                while n1 == n2:
                    n2 = random.randint(0, L-1)
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

def performStep2(state, L, alpha):

    for i in range(L):

        # Feedback effect
        r_feedback = random.random()
        if r_feedback < alpha:
            
            # Get random nucleosome
            n1 = random.randint(0, L-1)
            if state[n1] == 0:      # Will only change if unmethylated
                n2 = random.randint(0, L-1)
                while n2 == n1:
                    n2 = random.randint(0, L-1)

                if state[n2] != 0:  # 'n1' cannot change if 'n2' unmethylated
                    n3 = random.randint(0, L-1)
                    while n3 == n1 or n3 == n2:
                        n3 = random.randint(0, L-1)

                    if state[n2] == state[n3]:  # Both have same methylation -> change 
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
F = 10

tsteps = 1000000

# To store number of M/A and U at each timestep
nMU1 = np.empty((2,tsteps), dtype=int)
nMU2 = np.empty((2,tsteps), dtype=int)


alpha = fromFToAlpha(F)
state1 = getRandState(L)
state2 = getRandState(L)
for t in range(tsteps):
    state1 = performStep1(state1, L, alpha) # Model 1
    state2 = performStep2(state2, L, alpha) # Model 2

    nMU1[0,t] = sum(state1 == 0)    
    nMU1[1,t] = sum(state1 == 1)

    nMU2[0,t] = sum(state2 == 0)
    nMU2[1,t] = sum(state2 == 1)


np.savetxt("data/task2.1_F10.txt", nMU1, fmt="%d")
np.savetxt("data/task2.2_F10.txt", nMU2, fmt="%d")
