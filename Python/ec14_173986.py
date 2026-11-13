'''
Disciplina: EQ 220 - Métodos Numéricos
Atividade: Exercício Computacional 14 - Método de Runge-Kutta de 4ª Ordem
Nome:
    Daniel Mussato Campiotti
Data: 04/11/2026
Descrição: Slução de um p.v.i. de crescimento logístico utilizando-se do Método de Runge-Kutta de 4ª Ordem.
    dy / dt = 2 * y * (1 - y),
cuja solução analítica é dada por:
    y(t) = 1 / (1 + exp(-2 * t)).
comparando-se os resultado de Euler Explícito, EUler Implícito e Runge-Kutta de 4ª ordem.
    
'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# initial definitions:

def f(t, y):
    return 2 * y * (1 - y)

def runge_kutta_4(t0, yo, h, intervalo):
    t = np.arange(t0, intervalo + h, h)
    y = np.zeros(len(t))
    t[0] = t0
    y[0] = yo

    for i in range(len(t) - 1):
        k1 = h * f(t[i], y[i])
        k2 = h * f(t[i] + h / 2, y[i] + k1 / 2)
        k3 = h * f(t[i] + h / 2, y[i] + k2 / 2)
        k4 = h * f(t[i] + h, y[i] + k3)

        y[i + 1] = y[i] + (k1 + 2 * k2 + 2 * k3 + k4) / 6
        t[i + 1] = t[i] + h

    return t, y

def analytical_solution(t):
    return 1 / (1 + 4 * np.exp(-2 * t)) # where the initial condition is y(0) = 0.2

def euler_explícito(t0, y0, h, intervalo):
    t = np.arange(t0, intervalo + h, h)
    y = np.zeros(len(t))
    t[0] = t0
    y[0] = y0

    for i in range(len(t) - 1):
        y[i + 1] = y[i] + h * f(t[i], y[i])
        t[i + 1] = t[i] + h
    return t, y

def euler_implicito(t0, y0, h, intervalo):
    t = np.arange(t0, intervalo + h, h)
    y = np.zeros(len(t))
    y[0] = y0

    for i in range(len(t) - 1):
        erro = 1e0
        y_new = y[i]

        while (erro >= 1e-6):
            y_old = y_new
            y_new = y[i] + h * f(t[i + 1], y_old)
            erro = abs(y_new - y_old)

        y[i + 1] = y_new
        t[i + 1] = t[i] + h

    return t, y

def absolute_error(y_analytical, y_method):

    return np.linalg.norm(y_analytical - y_method)

# data for plotting:

t_rk, y_rk = runge_kutta_4(t0=0, yo=0.2, h=0.5, intervalo=8)
t_e, y_e = euler_explícito(t0=0, y0=0.2, h=0.5, intervalo=8)
t_i, y_i = euler_implicito(t0=0, y0=0.2, h=0.5, intervalo=8)
y_analytical = analytical_solution(t_rk)


# data dictionary with panda's DataFrame:

delta_rk = absolute_error(y_analytical, y_rk)
delta_e = absolute_error(y_analytical, y_e)
delta_i = absolute_error(y_analytical, y_i)


pd.set_option('display.float_format', lambda x: f'{x:.4e}')

data_normas = pd.DataFrame({
    'Método Numérico': ['Runge-Kutta (4th)', 'Euler Explícito', 'Euler Implícito'],
    '||erro||': [delta_rk, delta_e, delta_i]
})

# 4. Printando o resultado formatado
print("\n" + "_"*55)
print("     COMPARATIVO GLOBAL DE ERROS DOS MÉTODOS")
print("_"*55)
print(data_normas.to_string(index=False))
print("_"*55)

# graphic plotting:

t_analytical = np.linspace(0, 8, 100)
y_analytical = analytical_solution(t_analytical)

# plotting:
plt.figure(figsize=(10, 6))

plt.plot(t_analytical, y_analytical, color="black", linewidth=1.5, label="Solução Analítica")
plt.plot(t_rk, y_rk, marker="^", color="green", linestyle=":", markersize=8, label="Runge-Kutta (4th)")
plt.plot(t_e, y_e, ls=":", marker="o", color="purple", markersize=5, label="Euler Explícito")
plt.plot(t_i, y_i, ls=":", marker="s", color="red", markersize=5, label="Euler Implícito")

plt.xlabel("t")
plt.ylabel("y(t)")
plt.xlim(0, 8)
plt.ylim(0, 1.2)
plt.grid(True, linestyle=":", alpha=0.6)
plt.title("Comparação de Métodos Numéricos para Crescimento Logístico", fontsize=14, fontweight='bold')
plt.legend(fontsize=11)
plt.tight_layout()
plt.savefig("Graphics/ec14_173986_graphic.png", dpi=300, bbox_inches="tight")
plt.show()

