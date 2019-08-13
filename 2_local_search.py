# -*- coding: utf-8 -*-
'''
Created on Thu Jun 14 15:07:29 2018

@author: Maciej Workiewicz

The code has been tested on Python 2.7 and 3.6 and higher
'''

print('''
----------------------------------------------------
Running Module 2: Local search with long jumps
----------------------------------------------------
''')

import numpy as np
from os.path import expanduser  # new
import matplotlib.pyplot as plt

# Do not change these parameters--------------------------------------------
N = 6  #    set to N=6
i = 1000  # number of iterations, set to 1000 to save time
t = 50  # time periods set to 50 initially
# --------------------------------------------------------------------------

'''
You can experiment with different setting of the following two variables:
which_imatrix - select 1 for random IM; 2 modular; 3 nearly-mod; 4 diagonal

K - set to 2 for which_imatrix = 2, 3, and 4. For which_imatrix=1 you can choose
    K from 0 (no interactions) to N-1 (maximum interactions)
'''

# You can change those ---
which_imatrix = 1      # | type of the interaction matrix
K = 2                  # | number of interdependencies per decision variable
p_jump = 0.1           # | probability of a long jump in a given round
# ------------------------

if which_imatrix >1:  # to avoid a common mistake
    K = 2

# *** 1. LOAD THE NK LANDSCAPE FILE *****************************************

file_name = expanduser("~")
NK_landscape = np.load(file_name + '\\NK_workshop\\NK_land_type_' + str(which_imatrix) +
                       '_K_' + str(K) + '_i_' + str(i) + '.npy')

power_key = np.power(2, np.arange(N - 1, -1, -1))


# *** 2. LOCAL SEARCH WITH LONG JUMPS ***************************************

Output2 = np.zeros((i, t))

for i1 in np.arange(i):
    combination = np.random.binomial(1, 0.5, N)  # gen initial combination
    row = np.sum(combination*power_key)  # finding the address in the array
    fitness = NK_landscape[i1, row, 2*N]  # piggyback on work done previously
    max_fit = np.max(NK_landscape[i1, :, 2*N])
    min_fit = np.min(NK_landscape[i1, :, 2*N])
    fitness_norm = (fitness - min_fit)/(max_fit - min_fit)  # normalize 0 to 1
    for t1 in np.arange(t):  # time for local search
        Output2[i1, t1] = fitness_norm
        if np.random.rand() < p_jump:  # check whether we are doing a jump
            new_combination = np.random.binomial(1, 0.5, N)
        else:  # if not, then we simply search locally
            new_combination = combination.copy()
            choice_var = np.random.randint(N)
            new_combination[choice_var] = abs(new_combination[choice_var] - 1)
        row = np.sum(new_combination*power_key)
        new_fitness = NK_landscape[i1, row, 2*N]
        if new_fitness > fitness:  # if we have found a better combination
            combination = new_combination.copy()
            fitness = new_fitness.copy()
            fitness_norm = (fitness - min_fit)/(max_fit - min_fit)
        # otherwise all stays the same as in the previous round
Fitness2 = np.mean(Output2, axis=0)

# *** 3. PLOT ***************************************************************

# We will do a simple plot of the average fitness value over t periods
plt.figure(1, facecolor='white', figsize=(8, 6), dpi=150)  # for screens with
#          higher resolution change dpi to 150 or 200. For normal use 75.
plt.plot(Fitness2, color='green', linewidth=2, label='p_jump='+str(p_jump))
plt.ylim(0.5, 1)
plt.legend(loc=4,prop={'size':10})
plt.title('Results of local search', size=12)
plt.xlabel('time periods', size=12)
plt.ylabel('fitness', size=12)
print('Final fitness level for long jumps: ' + str(Fitness2[t-1]))

# END OF LINE
