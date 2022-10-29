import numpy as np
import random


def fromFToAlpha(F):
    return F/(1.+F)

def getRandState(L):
    state = np.empty(L, dtype=int)
    for i in range(L):
        state[i] = random.randint(-1,1)
    return state

def performStep(state, L, alpha):

    for i in range(L):

        # Feedback effect
        r_feedback = random.random()
        if r_feedback < alpha:
            n1 = random.randint(0, L-1)
            n2 = random.randint(0, L-1)
            while n1 == n2:
                n2 = random.randint(0, L-1)

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
F = [2, 4, 6]

tsteps = [10000, 10000, 1000000]

# To store number of M/A and U at each timestep
nMU = []
for i in range(len(F)):
    nMU.append(np.empty((2,tsteps[i]), dtype=int))

for idx, f in list(zip(range(len(F)), F)):
    alpha = fromFToAlpha(f)
    state = getRandState(L)
    for t in range(tsteps[idx]):
        state = performStep(state, L, alpha)
        nMU[idx][0,t] = sum(state == 0)
        nMU[idx][1,t] = sum(state == 1)


np.savetxt("data/task1_F2.txt", nMU[0], fmt="%d")
np.savetxt("data/task1_F4.txt", nMU[1], fmt="%d")
np.savetxt("data/task1_F6.txt", nMU[2], fmt="%d")
