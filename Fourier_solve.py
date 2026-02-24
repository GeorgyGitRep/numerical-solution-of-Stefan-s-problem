import numpy as np
import matplotlib.pyplot as plt




# Параметры задачи
L = 2# Длина стержня [м]
a2 = 1.2328042328042328e-06 #2.6*10e-6#4 Коэффициент температуропроводности [м^2/с]
T0 =-0.8#-Начальная температура 
T1 = -1.2   # Температура на левом конце 
T2 = 0.0    # Температура на правом конце 

# Параметры для ряда Фурье
def u_xt(x, t, N=50):
    w_x = T1 + (T2 - T1) * x / L

    v_xt = 0.0
    for n in range(1, N+1):
        C_n = (2/(np.pi*n)) * (T0 - T1 + ((-1)**n) * (T2 - T0))

        term = C_n * np.sin(np.pi * n * x / L) * np.exp(-a2 * (np.pi * n / L)**2 * t)
        v_xt += term
    
    return w_x + v_xt

x_points = np.linspace(0, L, 200)

times = [0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0,100,3600]

plt.figure(figsize=(12, 8))


for i, t in enumerate(times):
    u_values = [u_xt(x, t) for x in x_points]
    if i == len(times)-1:
        An_sol = u_values
        plt.plot(x_points, u_values, 
                 color='r', 
                 linewidth=2.5, 
                 label=f't = {t} с',
                 alpha=0.8)


plt.xlabel('x [м]', fontsize=12)
plt.ylabel('Температура, u(x,t) [°C]', fontsize=12)
plt.title('Распределение температуры в стержне, метод Фурье\n' + 
          f'Начальная температура: {T0}°C, Граничные условия: u(0,t)={T1}°C, u({L},t)={T2}°C', 
          fontsize=14)
plt.grid(True, alpha=0.3)
plt.legend(fontsize=11)


# Добавляем отметки начальной температуры
plt.axhline(y=T0, color='red', linestyle='--', alpha=0.3, label=f'Начальная: {T0}°C')

plt.tight_layout()
plt.show()



