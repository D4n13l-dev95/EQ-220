"""
Disciplina: EQ 220 - Métodos Numéricos
Atividade: Projeto Computacional 1
Nome:
    Daniel Mussato Campiotti
Data: 21/08/2026
Descrição: Implementar um código através de métodos numéricoss que encontre as soluções da equação de 
van der Walls para diferentes valores de T e P, para a água na interface líquido-vapor. 
O código deve plotar o gráfico de P x T. 
Além disso, o criitério de aceitação da solução do ponto (P,T) é baseado no critério dos valores de fugacidade (phi)
A princípio, o método numérico utilizado será Newton-Raphson.
"""

import numpy as np
import matplotlib.pyplot as plt

# definição de constantes
R = 8.314462e0       #cm3 · MPa · K−1 · mol−1
Tc = 647.1e0         # K
Pc = 22.06e0         # MPa

# Dados de referência do NIST WebBook para H2O (T [°C] -> Psat [kPa])
nist_data = {
    20: 2.339,  25: 3.169,  30: 4.247,  35: 5.629,  40: 7.385,
    45: 9.595,  50: 12.35,  55: 15.76,  60: 19.95,  65: 25.04,
    70: 31.20,  75: 38.58,  80: 47.41,  85: 57.87,  90: 70.18,
    95: 84.61, 100: 101.42
}

# definição de funções

# a e b para temp. crítica (eq. (2) do pc_1)
def calc_ab(Tc, Pc):
    a = 27 * R**2 * Tc**2 / (64 * Pc)
    b = R * Tc / (8 * Pc)
    return a, b

# A e B do eq. (6) depndentes de a e b 
def calc_AB(P, T, a, b):
    A = a * P / (R**2 * T**2)
    B = b * P / (R * T)
    return A, B

# eq. de van der Walls
def F(Z, A, B):

    Y = Z**3 - (1 + B) * Z**2 + A * Z - A * B

    return Y

# derivada da eq. de van der Walls
def dF(Z, A, B):

    dY = 3 * Z**2 - 2 * (1 + B) * Z + A

    return dY


# definação da fugacidade para qq. Z_liq ou Z_vap, basta camaá-los na função
def phi(Z, A, B):
    ln_phi = Z - 1 - np.log(Z - B) - A / Z
    return np.exp(ln_phi)


### listas ###

#lista de Temp no intervalo de 20 a 100 com incremento de 5

ls_T = np.arange(20, 105, 5) + 273.15    # lista de T de 20 a 100 com incremento de 5 (conversão para K)

ls_P = np.linspace(1.0e-4, 0.15, 1500)

ls_Z = np.linspace(1.0e-6, 1.5e0, 1000)
# N.B. - escolheu-se este intervalo, porque é um range que comporta a previsão experimental

# inicialização de a e b
valor_a, valor_b = calc_ab(Tc, Pc)

# print(valor_a, valor_b)
i = 0
# inicilaização de A e B

ls_raizes = []
for temperatura in ls_T:
    for pressao in ls_P:

        valor_A, valor_B = calc_AB(pressao, temperatura, valor_a, valor_b)

        raizes_desse_ponto = []
        for z in range(len(ls_Z) - 1):
            Z1 = ls_Z[z]
            Z2 = ls_Z[z + 1]
            F1 = F(Z1, valor_A, valor_B)
            F2 = F(Z2, valor_A, valor_B)
            if F1 * F2 < 0:
                raizes_desse_ponto.append((Z1, Z2))

        ls_raizes.append(raizes_desse_ponto)

print(len(ls_raizes))

        # Valor de Z 

