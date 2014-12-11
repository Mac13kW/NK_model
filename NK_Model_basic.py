print """ ----------------------------------------------------
NK model (Version 3.1)

----------------------------------------------------

by Maciej Workiewicz (2014)

----------------------------------------------------"""

# *** IMPORTS ***

import numpy as np
import itertools
from time import time

start = time()  # this will tell me how much time one run takes

# *** MODEL INPUTS ***

# SYSTEM INPUT
N = 4
i = 100  # number of landscapes to produce (use 100 or less for testing)

# USER INPUTS

print ('''
       Interaction matrix: 1 - random or 2 - custom
       ''')
which_matrix = int(raw_input("Choose interaction matrix (1 or 2): "))


# FUNCTIONS AND INTERACTION MATRIX


def matrix_rand(N, K):
    '''
    This function takes the number of elements and interdependencies and
    creates a random interaction matrix with a diagonal filled.
    '''
    Int_matrix_rand = np.zeros((N, N))
    for aa1 in np.arange(N):
        Indexes_1 = range(N)
        Indexes_1.remove(aa1)
        np.random.shuffle(Indexes_1)
        Indexes_1.append(aa1)
        Chosen_ones = Indexes_1[-(K+1):]  # this takes the last K+1 indexes
        for aa2 in Chosen_ones:
            Int_matrix_rand[aa1, aa2] = 1
    return(Int_matrix_rand)


if which_matrix == 1:
    K = int(raw_input("Input K (integer from 0 to N-1): "))

elif which_matrix == 2:  # custom
    K = 1  # remember to change that if you alter the pattern
    '''
    The diagonal of the interaction matrix is always 1
    '''
    Int_matrix = np.array([[1,1,0,0],\
                           [1,1,0,0],\
                           [0,0,1,1],\
                           [0,0,1,1]])

#ADDITIONAL FUNCTIONS


def powerkey(N):
    '''
    Used to find the location on the landscape for a given combination
    of the decision variables. Returns a vector with the powers of two
    (2^(N-1), ..., 1)
    '''
    Power_key = np.power(2, np.arange(N - 1, -1, -1))
    return(Power_key)


def nkland(N):
    '''
    Generates an NK landscape - an array of random numbers ~U(0, 1).
    '''
    NK_land = np.random.rand(2**N, N)
    return(NK_land)


def calc_fit(N, NK_land, inter_m, Current_position, Power_key):
    '''
    Takes landscape and a combination and returns a vector of fitness
    values (contribution value for each of for N decision variables)
    '''
    Fit_vector = np.zeros(N)
    for ad1 in np.arange(N):
        Fit_vector[ad1] = NK_land[np.sum(Current_position * inter_m[ad1]
                                         * Power_key), ad1]
    return(Fit_vector)


def comb_and_values(N, NK_land, Power_key, inter_m):
    """
    Calculates values for all combinations on the landscape.
    - the first N columns are for the combinations of N decision variables DV
    - the second N columns are for the contribution values of each DV
    - the next valuer is for the total fit (avg of N contributions)
    - the last one is to find out whether it is the local peak (0 or 1)
    """
    Comb_and_value = np.zeros((2**N, N*2+2))
    c1 = 0  # starting counter for location
    for c2 in itertools.product(range(2), repeat=N):
        '''
        this takes time so carefull
        '''
        Combination1 = np.array(c2)  # taking each combination
        fit_1 = calc_fit(N, NK_land, inter_m, Combination1, Power_key)
        Comb_and_value[c1, :N] = Combination1  # combination and values
        Comb_and_value[c1, N:2*N] = fit_1
        Comb_and_value[c1, 2*N] = np.mean(fit_1)
        c1 = c1 + 1
    for c3 in np.arange(2**N):  # now let's see if that is a local peak
        loc_p = 1  # assume it is
        for c4 in np.arange(N):  # check for the neighbourhood
            new_comb = Comb_and_value[c3, :N].copy()
            new_comb[c4] = abs(new_comb[c4] - 1)
            if ((Comb_and_value[c3, 2*N] <
                 Comb_and_value[np.sum(new_comb*Power_key), 2*N])):
                loc_p = 0  # if smaller than the neighbour then not peak
        Comb_and_value[c3, 2*N+1] = loc_p
    return(Comb_and_value)


# *** GENERAL VARIABLES AND OBJECTS ***

# For each iteration i

Power_key = powerkey(N)
NK = np.zeros((i, 2**N, N*2+2))

for i_1 in np.arange(i):
    '''
    First create the landscapes for the DDs and PDs
    '''
    NK_land = nkland(N)
    if which_matrix == 1:  # random
        Int_matrix = matrix_rand(N, K)

    NK[i_1] = comb_and_values(N, NK_land, Power_key, Int_matrix)

disc_loc = ('C:\\Output_folder\\')  # change that to match yours
np.save(disc_loc + 'NK_landscape_N_' + str(N) + '_K_' + str(K) +
        '_i_' + str(i) + '.npy', NK)
'''
This saves the landscape into a numpy binary file
'''

elapsed_time = time() - start
print(' time: ' + str("%.2f" % elapsed_time) + ' sec')

#END OF LINE
