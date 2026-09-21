"""
Disciplina: EQ 220 - Métodos Numéricos
Atividade: Exercício Computacional 7

Nome: Daniel Mussato Campiotti

Data: 17/09/2026

Descrição: Resolução de múltiplos balanços de massa
através do método iterativo de Jacobi. Implementação subsequente do método de Gauss-Seidel.
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


# tamanho do sistema

n = len(v)


# tolerância

tolerancia = 1e-7


# chute inicial

x_old = np.zeros(n)

x_old[0] = 10

# erro inicial

erro = 1


# método iterativo de Jacobi

m = 0  # counter

while erro >= tolerancia:

    x_new = np.zeros(n)

    i = 0

    while i < n:

        s = 0

        j = 0

        # termos antes da diagonal ("matriz L")

        while j < i:

            s = s + A[i][j] * x_old[j]

            j += 1

        # termos depois da diagonal

        j = i + 1

        while j < n:

            s = s + A[i][j] * x_old[j]

            j += 1

        # cálculo do novo x_i

        x_new[i] = (v[i] - s) / A[i][i]

        i += 1

    # cálculo do erro

    erro = np.sqrt(np.sum((x_new - x_old)**2))

    # atualização

    x_old = x_new.copy()

    m += 1


# apresentação dos resultados

print("\nSolução do sistema (Jacobi):")

for i in range(n):
    print(f"x{i + 1} = {x_old[i]:.6f}")

print(f"\nErro final = {erro:.6e}")

print(f"\nNúmero de iterações = {m}")




# método iterativo de Gauss-Seidel

# chute inicial

x_old = np.zeros(n)

x_old[0] = 10

erro = 1

m = 0  # counter

while erro >= tolerancia:

    x_new = np.zeros(n)

    i = 0

    while i < n:

        s = 0

        j = 0

        # termos antes da diagonal ("matriz L")

        while j < i:

            s = s + A[i][j] * x_new[j]
            j += 1

        # termos depois da diagonal

        j = i + 1

        while j < n:

            s = s + A[i][j] * x_old[j]

            j += 1

        # cálculo do novo x_i

        x_new[i] = (v[i] - s) / A[i][i]

        i += 1

    # cálculo do erro

    erro = np.sqrt(np.sum((x_new - x_old)**2))

    # atualização

    x_old = x_new.copy()

    m += 1


# apresentação dos resultados

print("\nSolução do sistema (Gauss-Seidel):")

for i in range(n):
    print(f"x{i + 1} = {x_old[i]:.6f}")

print(f"\nErro final = {erro:.6e}")

print(f"\nNúmero de iterações = {m}")

print("\nVê-se que o n° de iterações para o Método de Gauss-Seiel é inferior ao de jacobi:")
print("constrói-se a nova solução componente por componente e aproveita-se imediatamente aquilo que ")
print("se acabou de calcular.")