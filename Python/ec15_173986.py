'''
Disciplina: EQ 220 - Métodos Numéricos
Atividade: Exercício Computacional 15 - Problema de Valor de Contorno (p.v.c.) resolvido por 
Nome:
    Daniel Mussato Campiotti
Data: 04/11/2026
Descrição: Condução de calor em uma chapa bidimensional, com contorno:

    -superior => 100°C
    -inferior => 0°C
    -direito  => 50°C
    -esquerdo => 75°C
    
'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# algoritimo de contrução da matriz/sistema linear para solução do problema

n = 100

# um pouco distante daquele visto em sala, onde vimo a solução em 1D!
t = np.zeros((n,n))

# definição das bordas 
t[0,:] = 100.0 
t[-1, :] = 0.0
t[:, 0] = 75.0
t[:, -1] = 50.0

b = np.zeros(n * n) # T_i,j posições

i = 0
while (i < n):
    j = 0
    while (j < n):
        k = (i * n) + j
        b[k] = 4.0 * t[i, j]
        j += 1
    i += 1

a = np.eye(n * n) * 4.0

k = n + 1

while (k < (n - 1) * n):
    l = 1
    while (l < (n - 1) * n):
        if ((k == l) and (k % n != 0) and ((k + 1) % n != 0)):
            a[k, l - 1] = -1.0  #  esquerda
            a[k, l + 1] = -1.0  # direita
            a[k, l - n] = -1.0  # superior
            a[k, l + n] = -1.0  # baixo
        l += 1
    k += 1

T = np.linalg.solve(a, b)
T = T.reshape((n,n))

# plotting 
fig, ax = plt.subplots(figsize=(8, 8))

mapa_calor = ax.imshow(T[1:-1, 1:-1], cmap='coolwarm', vmin=0, vmax=100, interpolation='bilinear')

ax.set_xticks([])
ax.set_yticks([])

barra_lateral = plt.colorbar(mapa_calor, ax=ax, extend='both', shrink=0.9)
barra_lateral.set_label('Temperatura / °C', fontsize=14)

plt.tight_layout()
plt.savefig("Graphics/ec15_graphics.png")
plt.show()