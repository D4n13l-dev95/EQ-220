'''
Disciplina: EQ 220 - Métodos Numéricos
Atividade: Exercício Computacional 14 - Método de Runge-Kutta de 4ª Ordem
Nome:
    Daniel Mussato Campiotti
Data: 16/11/2026
Descrição: The main objective of this program is to elaborate a multiple solutions, for a system of o.d.e.'s.
By Runge-Kutta 4th order method, for a system of o.d.e.'s, with the following form:

    dx/dt = c1 * x
    dy/dt = (c2 - y) * c3 - c4 * y
    dz/dt = c5 * x - c6 * z 
    
'''

import numpy as np 
import matplotlib.pyplot as plt

# initial defiitions e.d.o.'s function: 

def system(t, state, c1, c2, c3, c4, c5, c6):

    x, y, z = state

    dx_dt = c1 * x 

    dy_dt = (c2 - y) * c3 - c4 * y

    dz_dt = c5 * x - c6 * z

    return np.array([dx_dt, dy_dt, dz_dt])


# Runge-Kutta 4th order method for a system of o.d.e.'s

def runge_kutta_step(f, t, state, h):

    k1 = f(t, state)
    k2 = f(t + 0.5 * h, state + 0.5 * h * k1)
    k3 = f(t + 0.5 * h, state + 0.5 * h * k2)
    k4 = f(t + h, state + h * k3)

    return state + h / 6 * (k1 + 2*k2 + 2*k3 + k4)

# bounary and initial conditions

t_start = 0e0
t_end = 10.0

# array for h variyng:

h_values = np.linspace(0.25, 1.0, num=4, endpoint=True)

def time_vector(t_start, t_end, h):
    time = np.arange(t_start, t_end + h/2, h)
    return time

# boundary conditions 

# x0 = 
# y0 = 
# z0 = 0e0

initial_condition = np.array([x0, y0, z0])

all_results = {}

# integrating on time 

for step in h:

    time_vec = time_vector(t_star, t_end, step)

    results = np.zeros((len(time_vec), 3))

    current_state = initial_conditions.copy()


    for i, t in enumerate(time_vec):
        results[i] = curent_state
        current_state = runge_kutta_step(system, t, current_state, h)

    all_results[step] = (time_vec, results)

plt.figure(figsize=(8, 5))
for step, (t_vec, res) in all_results.items():
    plt.plot(t_vec, res[:, 0], marker='o', markersize=3, label=f'h = {step}')

plt.title('Comparação da variável x(t) para diferentes passos h')
plt.xlabel('Tempo (t)')
plt.ylabel('x(t)')
plt.legend()
plt.grid(True)
plt.show()

