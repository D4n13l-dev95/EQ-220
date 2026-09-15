"""
Disciplina: EQ 220 - Métodos Numéricos
Atividade: Exercício Computacional 6

    Nome: Daniel Mussato Campiotti
    Nome: Maithe Andres Saad

Data: 09/09/2026

Descrição: Resolução de múltiplos balanços de massa através da decomposição LU.
"""

import numpy as np


# matriz A (coeficientes) e vetor v (termos independentes)

A = np.array([
    [-6.0, 0.0, 1.0, 0.0, 0.0],
    [3.0, -3.0, 0.0, 0.0, 0.0],
    [0.0, 1.0, -9.0, 0.0, 0.0],
    [0.0, 1.0, 8.0, -11.0, 2.0],
    [3.0, 1.0, 0.0, 0.0, -4.0]
])

v = np.array([-50.0, 0.0, -160.0, 0.0, 0.0])


# inicialização das matrizes L e U

n = len(v)

L = np.eye(n)
U = A.copy()


# eliminação de gauss para obter a matriz U
# os fatores de eliminação são armazenados na matriz L

h = 0
k = 0

while (h < n and k < n):

    i = h + 1

    while (i < n):

        # cálculo do fator de eliminação

        f = U[i][k] / U[h][k]

        # armazenamento do fator na matriz L

        L[i][k] = f

        # eliminação do elemento abaixo do pivô

        U[i][k] = 0.0

        j = k + 1

        while (j < n):

            U[i][j] = U[i][j] - U[h][j] * f

            j += 1

        i += 1

    h += 1
    k += 1


# substituição progressiva para resolver Lu = v

u = [0.0] * n

i = 0

while (i < n):

    s = 0.0

    j = 0

    while (j < i):

        s = s + L[i][j] * u[j]

        j += 1

    u[i] = (v[i] - s) / L[i][i]

    i += 1


# substituição reversa para resolver Ux = u

x = [0.0] * n

x[n - 1] = u[n - 1] / U[n - 1][n - 1]

i = n - 2

while (i >= 0):

    s = 0.0

    j = i + 1

    while (j < n):

        s = s + U[i][j] * x[j]

        j += 1

    x[i] = (u[i] - s) / U[i][i]

    i -= 1


# apresentação dos resultados

print("\nMatriz Triangular Inferior (L):")

for linha in L:
    print([round(valor, 4) for valor in linha])


print("\nMatriz Triangular Superior (U):")

for linha in U:
    print([round(valor, 4) for valor in linha])


print("\nVetor v:", [round(valor, 4) for valor in v])

print("\nVetor u:", [round(valor, 4) for valor in u])

# apresentação das concentrações obtidas

print("\nSolução do Sistema (x):", [round(valor, 4) for valor in x], ", i.e.:")

print(f"c1 = {x[0]:.4f} L/min")
print(f"c2 = {x[1]:.4f} L/min")
print(f"c3 = {x[2]:.4f} L/min")
print(f"c4 = {x[3]:.4f} L/min")
print(f"c5 = {x[4]:.4f} L/min")