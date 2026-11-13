"""
Disciplina: EQ 220 - Métodos Numéricos
Atividade: Projeto Computacional 1
Nome:
    Daniel Mussato Campiotti
Data: 21/08/2026
Descrição: Descobrir a raiz de uma função pelo método da Regula Falsi e da Bissecção

"""

import numpy as np

# Definição inicial da função

def f(x):
  epsilon = 2.0e-4
  D = 0.2e0
  Re = 1.2e5

  y = (np.sqrt(x)) * (-2 * np.log10(epsilon / (3.7e0 * D) + 2.51 / (Re * np.sqrt(x)))) - 1

  return  y

tolerancia = 1.0e-6   # definição da tolerância

# Inicialização de variáveis

a = 0.01e0
fa = f(a)

b = 0.5e0
fb = f(b)

if (fa * fb >= 0):

  print("Não há raiz no intervalo [", a, ",", b, "]")

  exit("Não há raiz no intervalo especificado.")

erro = abs(a - b)

c_anterior = b  # inicialização de c_anterior (i.e. próprio b inicial)

i = 0

# loop de iteração (Regula Falsi)
while (erro >= tolerancia):
  c = b - (f(b) * (b - a)) / (f(b) - f(a))
  fc = f(c)

  if (fa * fc >= 0):
    a = c
    fa = fc

  else:
    b = c
    fb = fc
   
  erro = abs(b - c_anterior)
  c_anterior = c
  i += 1

# Regula Falsi 
print(10*'---')
print("Convergiu por regula falsi!")
print("Raiz =", c)
print("Prova Real f(z) =", f(c))
print("Erro =", erro)
print("Número de iterações =", i)


# Bissecção

x = 0.01e0
fx = f(x)

y = 0.5e0
fy = f(y)

if (fx * fy >= 0):

  print("Não há raiz no intervalo [", x, ",", y, "]")

  exit("Não há raiz no intervalo especificado.")

erro = abs(x - y)

z_anterior = y  # inicialização de z_anterior (i.e. próprio y inicial)

j = 0

# loop de iteração (Bissecção)
while (erro >= tolerancia):

    z = (x + y) / 2.0e0
    
    fz = f(z)
    
    if (fx * fy >= 0):
        x = z
        fx = fz
    
    
    else:
        y = c
        fy = fz
     
    erro = abs(y - z_anterior)
    z_anterior = z
    j += 1

# Bissecção
print(10*'---')
print("Convergiu por bissecção!")
print("Raiz =", z)
print("Prova Real f(z) =", f(z))
print("Erro =", erro)   
print("Número de iterações =", j) 
print(10*'---')


# c. Comparação dos métodos
#  O número de iterações para o método da bissecção é maior do que o número de iterações para o método da regula falsi. 
# Isso ocorre porque o método da bissecção divide o intervalo ao meio a cada iteração, 
# enquanto o método da regula falsi utiliza uma abordagem mais eficiente, 
# ajustando o ponto de teste com base nos valores da função nos extremos do intervalo. 
# Portanto, a regula falsi tende a convergir mais rapidamente para a raiz em comparação com a bissecção.