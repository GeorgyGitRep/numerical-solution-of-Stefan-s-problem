import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.special import erf, erfc

def bisection_method(f, a, b, tol=1e-10, max_iter=100):
    fa = f(a)
    fb = f(b)
    
    if fa * fb > 0:
        raise ValueError(f"На границах интервала f(a)*f(b) > 0. f({a}) = {fa}, f({b}) = {fb}")
        iterations = 0
    for i in range(max_iter):
        c = (a + b) / 2
        fc = f(c)
        error = (b - a) / 2

        if error < tol or abs(fc) < tol:
            return c, i + 1

        if fa * fc < 0:
            b = c
            fb = fc
        else:
            a = c
            fa = fc
        
        iterations = i + 1
    c = (a + b) / 2
    return c, iterations

def equation(alpha):
    c = 2100 
    k1,  T1 = 2.330, -1.2 
    lmbda, rho = 334 * 1000, 900 
    a1 = (k1/(rho*c))
    beta = alpha / (2*a1)       
    left = 1/np.sqrt(np.pi) * (np.exp(-beta * beta)/erf(beta))   
    D = (lmbda * rho * a1*a1) / (k1 * T1)
    right = - D * beta
    return left - right 
alfa = bisection_method(equation, 0.00000001, 2, tol=1e-10, max_iter=300)[0]

Tleft = -1.2
c = 2100
k1,  T1 = 2.330, -0.8 
lmbda, rho = 334 * 1000, 900 
a1 = (k1/(rho*c))
B1 = -Tleft/(erf(alfa/(2*a1)))
t_end =  3600
X_end = alfa *np.sqrt(t_end)


X_end2 = np.sqrt(-2*t_end*a1*T1)
print(X_end)
print(X_end2)











