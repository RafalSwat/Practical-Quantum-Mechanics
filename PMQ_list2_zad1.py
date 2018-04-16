import matplotlib.pyplot as plt

# needs explicit subpackage import
from scipy.integrate import odeint
from scipy.optimize import fsolve, brentq
import numpy as np
import RKMethods

w = 4.0
a = -w*0.5
b = w*0.5

N = 500
z = np.linspace(a, b, N)

psi_0 = 0.
psit_0 = 1.

#search range
E = np.linspace(0, 50, 1000)
#potential
V = np.zeros(N)

#effective mass
m0 = 0.067
m = m0 * np.ones(N)


def calc_psi(z, psi, psit, V, m, E):
    dz = z[1]-z[0]
    for i in range(1, len(z)):
        psi[i] = psi[i-1] + dz*m[i]*psit[i-1]
        psit[i] = psit[i-1] + dz * 2.0 * (V[i] - E) * psi[i-1]
    return psi


psi = np.zeros(N)
psit = np.ones(N)

e_prev = 0.0
de_eps = 0.1
solutions = []
# Dla zadanego z góry spektrum energii
# sprawdzamy czy ktoras wartosc jest taka ze da 0 na drugim brzegu.
for e in E:
    p = calc_psi(z, psi, psit, V, m, e)

    previous_psi = p[-2]
    last_psi = p[-1]
    '''
    Filtrowanie danych.
    '''
    # Po pierwsze sprawdzamy czy zmienia się znak dwóch ostatnich wartości
    # Jesli tak to znaczy że "przeszło" przez 0
    if np.sign(previous_psi) + np.sign(last_psi) == 0:
        # Po drugie nie chcemy aby notowac kilkukrotnie te samy wartości
        de = e - e_prev
        if de > de_eps:
            print(previous_psi, last_psi, e, de)
            #kopia bo potem zmieniamy to co jest pod p (a to numpy ;))
            solutions.append( (e, p.copy()))
        e_prev = e


for s in solutions:
    _psi = s[1]
    dz = z[1] - z[0]
    dpsi = np.gradient(_psi, dz)
    e0 = s[0]
    for i in range(1000):
        e0 = e0 - _psi[-1]/dpsi[-1]
        print (e0)

for s in solutions:
    plt.plot(z, s[1], label=r'E=%0.2f' % s[0])

plt.grid()
plt.title("Wave functions for infinite potential well")
plt.legend()
plt.show()
