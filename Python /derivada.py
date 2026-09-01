"""
Exemple of derivative usage by bib. sympy 
"""

import sympy as sp
import numpy as np

x = sp.symbols('x')

# Define a função usando 'def'

def f(x):

    y = x**3 + sp.cos(x)

    return y

def df(x):  # Outra forma de definir de modo que se tem doois def e se chama cada uma separadamente

    dy = sp.diff(f(x), x)

    return dy



# Calcula a derivada normalmente
derivada = sp.diff(f(x), x)
print(derivada)  # Resultado: 3*x**2 - sin(x)

dy = df(x)
print(df(x))
 # Avalia a derivada em x=2