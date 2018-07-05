# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 15:07:29 2018

@author: Maciej Workiewicz
"""

print("""
----------------------------------------------------
Running Module 3: Decentralized local search
----------------------------------------------------
""")

import numpy as np
from os.path import expanduser  # new
import matplotlib.pyplot as plt

# Do not change these parameters----
N = 6  #                           |
i = 1000  # number of iterations   |
t = 50  # time periods             |
# ----------------------------------

'''
You can experiment with different setting of the following two variables:
which_imatrix - select 1 for random IM; 2 modular; 3 nearly-mod; 4 diagonal

K - set to 2 for which_imatrix = 2, 3, and 4. For which_imatrix=1 you can choose
    K from 0 (no interactions) to N-1 (maximum interactions)
'''

# You can change those
which_imatrix = 1  # | type of the interaction matrix
K = 2              # | number of interdependencies per decision variable
reorg = 50         # | reorganization: in which round we merge the two units
# --------------------

# *** 1. LOAD THE NK LANDSCAPE FILE *****************************************

file_name = expanduser("~")
NK_landscape = np.load(file_name + '\\NK_land_type_' + str(which_imatrix) +
                       '_K_' + str(K) + '_i_' + str(i) + '.npy')
# remember to change \\ into / on a Mac

power_key = np.power(2, np.arange(N - 1, -1, -1))


# *** 2. DECENTRALIZED LOCAL SEARCH *****************************************
'''
This simulation is based on the module 2.
'''
Output3 = np.zeros((i, t))

for i1 in np.arange(i):
    combination = np.random.binomial(1, 0.5, N)  # initial combination
    row = np.sum(combination*power_key)
    fitness = NK_landscape[i1, row, 2*N]
    max_fit = np.max(NK_landscape[i1, :, 2*N])  # use for normalization of perf
    min_fit = np.min(NK_landscape[i1, :, 2*N])  # ditto
    fitness_norm = (fitness - min_fit)/(max_fit - min_fit)
    for t1 in np.arange(t):
        Output3[i1, t1] = fitness_norm
        row = np.sum(combination*power_key)
        fitA = np.mean(NK_landscape[i1, row, N:int(N+N/2)])
        fitB = np.mean(NK_landscape[i1, row, int(N+N/2):int(N*2)])
        if t1 < reorg:  # for decentralized search
            '''
            The basic idea here is that we will split the vector of decision
            variables in two halfs {N0, N1, N2} and {N3, N4, N5}. We will then
            perform parallel search. We will locally search and compare
            resulting fitness changes separatelly.
            '''
            new_combination = combination.copy()
            new_combA = combination[:int(N/2)].copy()
            new_combB = combination[int(N/2):].copy()
            choice_varA = int(np.random.randint(0, int(N/2)))
            choice_varB = int(np.random.randint(0, int(N/2)))
            new_combA[choice_varA] = abs(new_combA[choice_varA] - 1)
            new_combB[choice_varB] = abs(new_combB[choice_varB] - 1)
            new_combination[:int(N/2)] = new_combA.copy()
            new_combination[int(N/2):] = new_combB.copy()
            
            row = np.sum(new_combination*power_key)  # find address for new comb
            new_fitA = np.mean(NK_landscape[i1, row, N:(int(N+N/2))])  # fitness goal 1
            new_fitB = np.mean(NK_landscape[i1, row, (int(N+N/2)):int(N*2)])  # fitness goal 2
            if new_fitA > fitA:
                combination[:int(N/2)] = new_combA.copy()
            if new_fitB > fitB:
                combination[int(N/2):] = new_combB.copy()
            row = int(np.sum(combination*power_key))
            fitness = np.mean(NK_landscape[i1, row, N:2*N])  # final fitness
        else:  # otherwise we centralize the local search
            new_combination = combination.copy()
            choice_var = np.random.randint(N)
            new_combination[choice_var] = abs(new_combination[choice_var] - 1)
            row = np.sum(new_combination*power_key)
            new_fitness = NK_landscape[i1, row, 2*N]
            if new_fitness > fitness:
                combination = new_combination.copy()
                fitness = new_fitness.copy()
        fitness_norm = (fitness - min_fit)/(max_fit - min_fit)
        
Fitness3 = np.mean(Output3, axis=0)

# *** 3. PLOT ***************************************************************

plt.figure(1, facecolor='white', figsize=(8, 6))
plt.plot(Fitness3, color='blue', linewidth=2, label='centr in t: ' + str(reorg))
plt.ylim(0.5, 1)
plt.legend(loc=4,prop={'size':10})
plt.title('Results of local search', size=12)
plt.xlabel('time periods', size=12)
plt.ylabel('fitness', size=12)

print('For IM=' + str(which_imatrix) + ' and reorg in round: ' + str(reorg))
print('Final fitness level for decentr/centr: ' + str(Fitness3[t-1]))



