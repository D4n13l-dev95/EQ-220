"""
Disciplina: EQ 220 - Métodos Numéricos
Atividade: Projeto Computacional 1
Nome:
    Daniel Mussato Campiotti
Data: 28/08/2026
Descrição: Interpolação polinomial para a modelagem do calor específico (à pressão contante) Cp do 
do sólido de Ni. 

    T (K) = 1500, 1900, 2300, 2700

    Cp (j.mol⁻¹.K⁻¹) = 29.4756, 32.160, 35.6409, 42.2840

Spline quadrático aos dados da EC9.
    
"""

import numpy as np
import matplotlib.pyplot as plt

X = np.array([1500, 1900, 2300, 2700], dtype=float) # values of T (K)

Y = np.array([29.4756, 32.160, 35.6409, 42.2840], dtype=float) # values of Cp (j.mol⁻¹.K⁻¹)

n = len(X) # vector's lenght

b = [] # coefficientes of the polynomial interpolation (go to: line 107)
