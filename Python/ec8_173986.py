"""

Disciplina: EQ 220 - Métodos Numéricos
Atividade: Projeto Computacional 1
Nome:
    Daniel Mussato Campiotti
Data: 28/08/2026
Descrição: Solução de sistema não linear, dado por:

| x*x + y*y + z*z = 6
| x + y - z = 1
| x * y * z = 1 

pelo Método de Newton Multidimensional, em que a solução é dada por 

    x_(k+1) = x_k - J⁻¹ * f(x_k)

"""

import numpy as np

# def. da função f isolada 
def f(v):
    x, y, z = v
    return np.array([
        x**2 + y**2 + z**2 - 6.0,
        x + y - z - 1.0,
        x * y * z - 1.0
    ])

# definição da matriz Jacobiana (J) simples
def J(v):
    x, y, z = v
    return np.array([
        [2.0 * x, 2.0 * y, 2.0 * z],
        [1.0, 1.0, -1.0],
        [y * z, x * z, x * y]
    ])

# parâmetros pré-iterativos   
x_old = np.array([2.0, 1.0, 1.0])

tolerance = 1.0e-7

norm = np.linalg.norm(x_old) # norma do vetor como parâmetro de erro 

n = 0 # counter

while (norm >= tolerance):

    J_k = J(x_old) # cal. inicial do jacobiano 

    J_inv = np.linalg.inv(J_k) # J -> J⁻¹

    x_new = x_old - J_inv @ f(x_old) # i.e. x_(k+1) = x_k - J⁻¹ * f(x_k)

    fx = f(x_new)

    norm = np.linalg.norm(fx) 
    # quanto mais próximo de 0 menor é o "ruído" do vaor obtido. Daí o utilisar da norma(f) como parâmetro de erro.

    x_old = x_new

    n += 1

print("Solução encontrada, pelo Método de Newton Multidimensional, em", n, "iterações:")
print(f"x = {x_new[0]:3f}")
print(f"y = {x_new[1]:3f}")
print(f"z = {x_new[2]:3f}")