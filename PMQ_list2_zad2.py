import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.optimize import newton

'''defining virables'''

m0 = 9.10938356
mw = 0.067 * m0
mb = 0.10 * m0
E = 0.08
V = 0.334
hbar = 1.0545718
width = [3,5,7,9]
b = 4

#g(alpha) = 2pi*alpha/(b+w)'''
#-1<alpha<1'''

'''defining wave func'''

def f(E, alpha):
    k = np.sqrt(2 * mw * E / hbar**2)
    v = np.sqrt(2 * mb * (V - E) / hbar**2)
    c = (mb**2 * k**2 - mw**2 * v**2) / (2 * mb * mw * k * v)
    A = np.cos(w * k) * np.cosh(b * v)
    B = np.sin(w * k) * np.sinh(b * v)
    return A - c * B - np.cos(alpha * 2 * np.pi)

'''newton for alpha(-1:1)'''
#for diff. width

af = np.linspace(-1,1, 200)
for w in width:
    e = []
    for alpha in af:
        root = newton(f, 0.08, args = (alpha,))
        e.append(root)
    plt.plot(af, e, label = w)
plt.grid()
plt.xlabel("wave vector g")
plt.ylabel("energy E")
plt.legend()    
plt.show()


