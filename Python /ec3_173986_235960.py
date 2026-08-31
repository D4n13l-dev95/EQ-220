'''
Disciplina: Métodos Numéricos - EQ220

Atividade: Exercício Computacional 3

Nomes:
    - Daniel Mussato Campiotti (173986)
    - Rafael Mazoli Kaysel Peres (235960)

Data: 24/08/2026

Datas de modificação : 26/08/2026 
(Implementação do item b. - Gráfico de [A] em função do parâmetro F;
Implementação de comentários e melhorias no código)

Descrição: O código a seguir resolve um balanço de massa para um reator (CSTR) em estado estacionário
 para a espécie A por Regula Falsi. O que se quer é saber a concentração final de A em função do 
 parâmetro F. O código também plota a concentração final de A ([A]) em função do parâmetro F.

'''

import numpy as np
import matplotlib.pyplot as plt

# definição inicial da função
def f(C_a, F, kappa = 4.0e0, K = 1.0e0, C_in = 2.0e0, V = 0.5e0): 
  # parece haver algum ganho de desempenho em definir os parâmetros constants como argumentos da função.

  y = F * (C_in - C_a) - V * (kappa * C_a / (1 + K * (C_a * C_a)))

  return y

tolerancia = 1.0e-6   # definição da tolerância

# inicialização de variáveis

x = 0.0e0           # mínimo possível em se tratando de concentração (i.e. C_a ~ 0, num ideal)
fx = f(x , F = 0.2e0)

y = 2.0e0           # máximo possível em se tratando de concentração (i.e. C_a < C_in)
fy = f(y , F = 0.2e0)

if (fx * fy >= 0):

  exit("Não há raiz no intervalo [", x, ",", y, "]")


erro = abs(x - y)

z_anterior = y  # inicialização de z_anterior (i.e. próprio y inicial)

i = 0

# loop de iteração (Regula Falsi para F = 0.2)
while (erro >= tolerancia):
  z = y - (f(y, F = 0.2e0) * (y - x)) / (f(y, F = 0.2e0) - f(x, F = 0.2e0))

  fz = f(z, F = 0.2e0)

  if (fx * fz >= 0):
    x = z
    fx = fz


  else:
    y = z
    fy = fz
   
  erro = abs(z - z_anterior)

  z_anterior = z

  i += 1

# Regula Falsi 
print()
print(5*'~~',"Resultados",5*'~~')
print("Convergiu por regula falsi!")
print("Raiz =", z)
print("Erro =", erro)
print("Número de iterações =", i)
print(11*'~~~')

# item b. (O  gráfico de [A] em função do parâmetro F)
F_valores = np.linspace(0.2, 2.0, 100)
Ca_valores = []

# laço de iteração para cada valor de F (i.e. novo Regula Falsi para cada valor de F)
for F in F_valores:
  x = 0.0e0               # mínimo possível em se tratando de concentração (i.e. C_a ~ 0, num ideal)
  fx = f(x , F = F)

  y = 2.0e0               # máximo possível em se tratando de concentração (i.e. C_a < C_in)
  fy = f(y , F = F)

  if (fx * fy >= 0):
    print("Não há raiz no intervalo [", x, ",", y, "]")
    exit("Não há raiz no intervalo especificado.")

  erro = abs(x - y)

  z_anterior = y  # inicialização de z_anterior (i.e. próprio y inicial)

  i = 0

  while (erro >= tolerancia):
    z = y - (f(y, F) * (y - x)) / (f(y, F) - f(x, F))
    fz = f(z, F)

    if (fx * fz >= 0):
      x = z
      fx = fz
    else:
      y = z
      fy = fz

    erro = abs(z - z_anterior)
    z_anterior = z
    i += 1

  Ca_valores.append(z)

plt.figure(figsize=(9, 5))

plt.plot(
  F_valores, 
  Ca_valores, 
  label = 'Concentração $[A]$ em função de $F$', 
  color = 'blue'
  )

plt.title(
  'Concentração $[A]$ em Função da Vazão Volumétrica $F$ (Regula Falsi)',
  fontsize=12
  )

plt.xlabel(
  'Vazão Volumétrica $F$ $(m^3 \cdot h^{-1})$',
  fontsize = 10
  )

plt.ylabel(
  'Concentração $[A]$ $(mol \cdot m^{-3})$', 
  fontsize = 10
  )
# Formatação texto gráfico, constantes do problema.

k = 4.0
kazao = 1.0
C_inicial = 2.0
Vol = 0.5

texto_grafico = (
  "Constantes do problema\n"
  f"$\kappa =${k}$\;h^{{-1}}$\n" 
  f"$K = ${kazao}$\;m^{6} \cdot mol^{{-2}}$\n"
  f"$C_i =$ {C_inicial} $\;mol \cdot m^{{-3}}$\n" 
  f"$V = ${Vol}$\;m^3$"
  )

plt.text(1.5, 0.4, texto_grafico, fontsize = 10, bbox = dict(facecolor = 'white', alpha = 0.5, edgecolor = 'black'))

plt.grid(
  True, 
  linestyle = ':', 
  alpha = 0.7
  )

plt.legend(
  loc='upper left', 
  fontsize=10
  )

plt.show()