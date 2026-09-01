"""
Disciplina: EQ 220 - Métodos Numéricos
Atividade: Projeto Computacional 1
Nome:
    Daniel Mussato Campiotti
Data: 01/09/2026
Descrição: Discover the root of a function using the Newton-Raphson method.
After that, discover the roots in the interval from -1 to 5 
and, finally, plot the graphics of the function and its derivative.

"""

import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

# Definição inicial da função

def f(x, a = 5, b = 6, c = 1):

  y = x**3 - a * (x*x) + b * x - c  
  
  return  y


def df(x, a = 5, b = 6):
  dy = 3 * (x**2) - 2 * a * x + b
  return dy  

tolerance = 1.0e-7   # definição da tolerância
i = 0                # contador de iterações

# initialization of the first guess

x_new = 0.0e0       
y = f(x_new)
dy = df(x_new)
erro = 1.0e0

while erro > tolerance:
    x_new = x_new - y / dy
    y = f(x_new)
    dy = df(x_new)
    x_old = x_new
    erro = abs(x_new - x_old)
    i += 1

print("Root finded: ", x_new)

# interval from -1 to 5 with step of 0.1

interval = np.arange(-1, 5, 0.1)
j = 0  
ls_root = []  # list to store the roots found in the interval

# second iteration to calculate the roots in the interval

for j in range(len(interval) - 1): # evaluating intervals with roots
    x1 = interval[j]

    x2 = interval[j + 1] 

    f1 = f(x1)

    f2 = f(x2)

    if f1 * f2 < 0:  # check if the function changes sign in the interval
        print("Root found in the interval: [", x1, ",", x2, "]")

        ls_root.append(x1)

        ls_root.append(x2)
    else:           # exit the loop if no root is found in the interval
      continue 

for k in range(0, len(ls_root), 1):  # finding the roots in the intervals found
    x_new = ls_root[k]

    y = f(x_new)

    dy = df(x_new)

    erro = 1.0e0

    while erro > tolerance:
        x_new = x_new - y/dy

        y = f(x_new)

        dy = df(x_new)

        x_old = x_new

        erro = abs(x_new - x_old)

        i += 1

    print("Root finded: ", x_new)
    print("Number of iterations: ", i)

# plotting the function and its derivative if you like to!

# y_valores = [f(x) for x in interval]
# dy_valores = [df(x) for x in interval]

# plt.figure(figsize=(10, 6))

# plt.plot(interval, 
#          y_valores, 
#          label='f(x)'
#          )

# plt.plot(interval, 
#          dy_valores, 
#          label="$df/dx$"
#          )

# plt.axhline(0,
#             color= 'black', 
#             linewidth=0.5
#             )

# plt.axvline(0, 
#             color='black', 
#             linewidth=0.5
#             )


# plt.xlabel('$x$')
# plt.ylabel('$y$')

# plt.title('$Function\;and\; Derivative$')

# plt.grid(True)

# plt.legend()

# plt.show()
