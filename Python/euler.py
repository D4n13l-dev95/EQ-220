'''
Disciplina: EQ 220 - Métodos Numéricos
Atividade: Exercício Computacional 12
Nome:
    Daniel Mussato Campiotti
Data: 27/10/2026
Descrição: Resolução do seguinte problema da valor inicial:
    dy / dt = -y³
em que,
    y(0) = 1.
    Pelo Método de Euler Explícito.
'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def euler_explicito(t_0, y_0, h):
    lista_t = [t_0]
    lista_y = [y_0]

    t_old = t_0
    y_old = y_0
    
    for i in range(100):
        y_new = y_old - k * h * y_old
        t_new = t_old + h
        lista_t.append(t_new)
        lista_y.append(y_new)
        t_old = t_new
        y_old = y_new

    return lista_t, lista_y



def euler_implicito_v1(t_0, y_0, h):
    lista_t = [t_0]
    lista_y = [y_0]

    t_old_e = t_0
    y_old_e = y_0
    for i in range(100):

        erro_nr = 1
        y_old_nr = y_old_e
        while erro_nr >= 1e-6:
            g = y_old_nr + h * k * y_old_nr - y_old_e
            dg = 1.0 + h * k

            y_new_nr = y_old_nr - g / dg

            erro_nr = abs(y_new_nr - y_old_nr)

            y_old_nr = y_new_nr

        y_new_e = y_new_nr
        t_new_e = t_old_e + h

        lista_t.append(t_new_e)
        lista_y.append(y_new_e)

        t_old_e = t_new_e
        y_old_e = y_new_e

    return lista_t, lista_y




t_0 = 0.0
y_0 = 1.0

lista_h = [0.2, 0.4, 0.6, 0.8]
k = 2.0

fig, axs = plt.subplots(2, 2)

for h, ax in zip(lista_h, axs.flatten()):
    t_e, y_e = euler_explicito(t_0, y_0, h)
    t_i, y_i = euler_implicito_v1(t_0, y_0, h)

    ax.plot(t_e, y_e, ls="--", marker="o", label="Euler explícito")
    ax.plot(t_i, y_i, ls=":", marker="s", label="Euler implícito")
    ax.set_ylabel("y")
    ax.set_xlabel("x")
    ax.set_title(f"h = {h}")
    ax.legend()
    ax.set_xlim(0, 5)
plt.show()