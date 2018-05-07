import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


'''
-------------------------------------------------------------------------------
Define Exercise 1
-------------------------------------------------------------------------------
'''

def exercise_1():

# declaration of variables
    
    me = 9.109383e-31                 # <- electron mass [kg]
    mx = me * 0.067                   # <- effective mass [kg]
    hbar = 1.053571e-34               # <- planck constant devide by 2pi [Js]
    eV = 1.602176e-19                 # <- electronovolt [J]
    n = 4                             # <- numbers of energy levels
    
# Define a loop over first 4 energy level, for well width range (-14,15)
    
    E = np.empty([n, 29])             # <- array thats holding the values of energy
    for ni in range(1, n + 1):        # <- iterating over energy lvls
        for wi in range(-14, 15):     # <- iterating in range (1,30)
            if wi == 0:               # <- we canot devide by 0,
                Ei = 100              #    so when we are in the middle, 
                E[ni-1][wi+15-1] = Ei #    we givethe value of energy by hand
            else:
                Ei = (hbar**2 * np.pi**2 * ni**2)/ (8 * mx * (wi/2)**2 * eV * 1e-9**2)  
                E[ni-1][wi+15-1] = Ei
    E = np.transpose(E)
    

# Ploting (using DataFrame from pandas)

    E=pd.DataFrame(E, columns=('n1','n2','n3','n4')) 
    plt.plot(range(-14, 15), E['n1'], 'ro')
    plt.plot(range(-14, 15), E['n2'], 'g+')
    plt.plot(range(-14, 15), E['n3'], 'y^')
    plt.plot(range(-14, 15), E['n4'], 'bs')
    plt.xlabel('Width (nm)')
    plt.ylabel('Energy (eV)')
    plt.xlim(-18, 18)
    plt.show()

'''
-------------------------------------------------------------------------------
Define Exercise 2
-------------------------------------------------------------------------------
'''

def exercise_2():

# declaration of variables
    
    me = 9.109383e-31            # <- electron mass [kg]
    hbar = 1.053571e-34          # <- planck constant devide by 2pi [Js]
    eV = 1.602176e-19            # <- electronovolt [J]
    V = [0.5, 1 , 2]             # <- potential levls
    E = np.zeros([8,3])          # <- array thats holding the values of energy
    j = 0                        # <- iterator over the array (incremented later in kode)
    for i in V:                  # <- iterating over potential lvls
        a = 1                    # <- a is just a sneaky way to break a loop :)
        ni = 1
        while a == 1:
            Ei = (hbar**2 * np.pi**2 * ni**2) / (2 * me * 0.067 * 10**2 * eV * 1e-9**2)
            if Ei < i:           # <- Enegry must by smaller then barrier V 
                E[ni-1][j] = Ei
                ni = ni + 1
            else:
                a = 0
        j = j + 1
    PLOT(E)
    
# define insane ploting function
    
def PLOT(E):
    
    plt.plot([0, 2], [E[0][0], E[0][0]], color='r', linestyle='--', linewidth=2)
    plt.plot([0, 2], [E[1][0], E[1][0]], color='r', linestyle='--', linewidth=2)
    plt.plot([2, 4], [E[0][1], E[0][1]], color='c', linestyle=':', linewidth=2)
    plt.plot([2, 4], [E[1][1], E[1][1]], color='c', linestyle=':', linewidth=2)
    plt.plot([2, 4], [E[2][1], E[2][1]], color='c', linestyle=':', linewidth=2)
    plt.plot([2, 4], [E[3][1], E[3][1]], color='c', linestyle=':', linewidth=2)
    plt.plot([4, 6], [E[0][2], E[0][2]], color='b', linestyle='-.', linewidth=2)
    plt.plot([4, 6], [E[1][2], E[1][2]], color='b', linestyle='-.', linewidth=2)
    plt.plot([4, 6], [E[2][2], E[2][2]], color='b', linestyle='-.', linewidth=2)
    plt.plot([4, 6], [E[3][2], E[3][2]], color='b', linestyle='-.', linewidth=2)
    plt.plot([4, 6], [E[4][2], E[4][2]], color='b', linestyle='-.', linewidth=2)
    plt.plot([0, 6], [0, 0], color='k', linestyle='-', linewidth=3)
    plt.ylabel('energy (eV)')
    plt.show()


'''
-------------------------------------------------------------------------------
Define Exercise 3
-------------------------------------------------------------------------------
'''

def exercise_3():
    
# declaration of variables
    
    m0 = 9.10938356            # <- electron mass
    mw = 0.067 * m0            # <- electron mass in the well
    mb = 0.10 * m0             # <- electron mass in the barier
    V0 = 0.334                 # <- potential energy of barier
    hbar = 1.0545718           # <- planck constant devide by 2pi 
    w = [3, 5, 8]              # <- different well width
    
    for wi in w:
        T = []                 # <- array thats holding the values of transition
        E = []                 # <- array thats holding the values of energy
        for e in np.linspace(0.01, 1.6, 10000):
            if e == V0:
                continue
            elif e < V0:
                E.append(e)
                k1 = np.sqrt((2 * mw * e) / (hbar**2))
                k2 = np.sqrt((2 * mb * (V0 - e)) / hbar**2)
                t = 1 / (1 + 0.25 * (k1 / k2 + k2 / k1)**2 * (np.sinh(k2 * wi))**2)
                T.append(t)
            else:
                E.append(e)
                k1 = np.sqrt((2 * mw * e) / (hbar**2))
                k2 = np.sqrt((2 * mb * (e - V0)) / hbar**2)
                t = 1 / (1 + 0.25 * (k1 / k2 - k2 / k1)**2 * (np.sin(k2 * wi))**2)
                T.append(t)
        plt.plot(E, T, label = "w =%s nm" %wi)
    plt.legend()
    plt.grid()
    plt.xlabel('E [eV]')
    plt.ylabel('Transmission coefficient \nT [0-1]')
    plt.show()

#exercise_1()
#exercise_2()   
#exercise_3()   
