# importing libraries
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
%matplotlib inline

# declaration of variables
me = 9.109383e-31   # kg  
hbar = 1.053571e-34 # Js
eV = 1.602176e-19   # J
mw = me * 0.067     # masa w studni
mb = me * 0.1       # masa w barierze
V = 0.334 * eV      # potencjal
E = 0.01 * eV       # energia czastki
b = 10 * 1e-9       # szerokosc bariery
w = 10 * 1e-9       # szerokosc studni
psi = []            # warosci f. falowej
ksi = []            # wartosci f.pomocniczej
iterator = []       # war. iteratora
psi.append(0)       # war. poczatkowe psi[0]=0
ksi.append(1)       # war. poczatkowe ksi[0]=0
iterator.append(0)  # war. poczatkowe i[0]=0

for z in range(0,1000) :
    if 300 < z < 700:
        psi.append(psi[z] + mb * ksi[z] * z)
        ksi.append(ksi[z] + (2/hbar) * (V - E) * psi[z] * z)
        iterator.append(z)
    else:
        psi.append(psi[z] + mw * ksi[z] * z)
        ksi.append(ksi[z] + (2/hbar) * (0 - E) * psi[z] * z)
        iterator.append(z)
    
plt.plot(iterator, psi)

