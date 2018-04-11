# deklaracja bibliotek
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import fsolve

# declaration of variables
me = 9.109383e-31   # kg  
hbar = 1.053571e-34 # Js
eV = 1.602176e-19   # J
h = 1e-12           # krok 
mw = me * 0.067     # masa w studni
mb = me * 0.1       # masa w barierze
V = 0.334 * eV      # potencjal
E0 = 0.1 * eV       # energia czastki
b = 10 * 1e-9       # szerokosc bariery
w = 10 * 1e-9       # szerokosc studni
fV = []             # wartosci bariery do wykresu
psi = []            # warosci f. falowej
ksi = []            # wartosci f.pomocniczej
iterator = []       # war. iteratora
psi.append(0)       # war. poczatkowe psi[0]=0
ksi.append(1)       # war. poczatkowe ksi[0]=0
iterator.append(0)  # war. poczatkowe i[0]=0
fV.append(V)        # war. poczatkowe fV[0]=0

# definiuje metode strzalow
def f_psi(E):
    h = 1e-12
    mw = me * 0.067 
    mb = me * 0.1
    hbar = 1.053571e-34
    
    for i in range(1,100):
        if 30 < i < 70:
            iterator.append(z * 10 * 10e-9)
            psi.append(psi[i-1] * h * mw * ksi[i-1])
            ksi.append(ksi[i-1] * h * (2.0/hbar) * (0 - E) * psi[i-1])
            fV.append(0)
        else:
            iterator.append(z * 10 * 10e-9)
            psi.append(psi[i-1] * h * mb * ksi[i-1])
            ksi.append(ksi[i-1] * h * (2.0/hbar) * (V - E) * psi[i-1])
            fV.append(V)
        return psi 
    
# definiuje pochodna psi
def df_psi(x):
    h = 1e-12 
    xh=x*h
    return 0.5*(f_psi(x+xh)-f_psi(x-xh))/xh

# definiuje metode Newtona
def Newton(x):
    eps = 10e-14
    x1 = x
    x0 = 0
    i = 0
    lmax = 500
    while np.fabs((x0 - x1)/x1) and (i < lmax):
        x0 = x1
        x1 = x0 - f_psi(x0)/df_psi(x0)
        i = i + 1
    if i < lmax:
        return x1
    else:
        return -100.0
    
# szukanie warosci wlasnych
max = 10e7
E0 = 0.01
for i in range(1, 334):
    E0 = E0 + i/1000
    Ew = Newton(E0)
    print(Ew)
