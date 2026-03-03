from matplotlib.pyplot import plot as plt
import numpy as np

Tleft = -1.2
Tright = 0.0
newLenght = []

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


N = 100
Tstart = ([-0.8]*N)
Tstart[0] = Tleft
Tstart[-1] = Tright
H = ([0.02]*(N-2))
    
#Tstart = [Tleft,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,Tright]
#H = [0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,] #0.02
print(sum(H))
lh = len(H)
tau = 20 #0.05
Tend = 3600
TN = int(Tend/tau)
h = 0.2
lamda =   1.2328042328042328e-06#  2.6*10e-6 #1.81976456306803*10e-6


T = Tstart.copy()
for j in range(TN):
    print(f'{j} из {TN}')
    alfa = [0]
    beta = [Tleft]
    #print(f"i:0 alfa = {alfa[0]} beta = {beta[0]}")
    
    for i in range(1,len(T)-2):
        a = 1/(H[i-1])
        b = 1/(H[i])
        c = a + b + (H[i-1] + H[i])/(2*lamda*tau)
        f = T[i] * ((H[i-1] + H[i])/(2*lamda*tau))
        beta.append((f+a*beta[-1])/(c-a*alfa[-1]))
        alfa.append((b)/(c-a*alfa[-1]))
        #print(f"i:{i} alfa = {alfa[i]} beta = {beta[i]}")
    
    
    def f(x):
        G = -334 * 1000
        a = 1/(H[-1])
        b = 1/(x)
        c = a + b + (H[-1] + x)/(2*lamda*tau)
        f = T[-2] * ((H[-1] + x)/(2*lamda*tau))
        Left = x*x * (c - a * alfa[-1])
        Right = tau/G*(f + a * beta[-1])
        return Left - Right
    
    #R = np.linspace(-0.01, 0.1,30)
    #U = [f(x) for x in R]
    #plt(R,U)
    x = bisection_method(f, 0.00000000000000001, 0.1, tol=1e-16, max_iter=300)[0] 
    #fsolve(f, 0)[0] #fsolve(f, 0)[0] 
    H.append(x)
    #print(f'H = {H}')
    newLenght.append(x) 
    
    a = 1/(H[i-1])
    b = 1/(H[i])
    c = a + b + (H[i-1] + H[i])/(2*lamda*tau)
    f = T[i] * ((H[i-1] + H[i])/(2*lamda*tau))
    beta.append((f+a*beta[-1])/(c-a*alfa[-1]))
    alfa.append((b)/(c-a*alfa[-1]))
    #print(f'j:{j}   newLenght = {newLenght}') 
    
    Tnew=[Tright]
    #print(len(T)," ",len(alfa))
    for k in range(len(T)-2,-1,-1):
        Tnew.append(alfa[k]*Tnew[-1] + beta[k])
    Tnew = Tnew[::-1]
    Tnew.append(Tright)
    T = Tnew.copy()
    
X=[0]
print(H)
sums = 0
for i in range(len(H)):
    sums = sums + H[i]
    print(sums)
    X.append(sums)
#plt(X,T[0:-1])
print(len(X))
sums = [0]
for i in newLenght:
    sums.append(sums[-1]+ i)
Time = np.linspace(0, Tend,TN)
print(sums[-1])
plt(Time,sums[:-1],label=f't = {Tend} с')
