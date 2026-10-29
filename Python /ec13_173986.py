'''
Disciplina: EQ 220 - Métodos Numéricos
Atividade: Exercício Computacional 13
Nome:
    Daniel Mussato Campiotti
Data: 28/10/2026
Descrição: 
    Resolução numérica de um Problema de Valor Inicial (PVI) definido por 
    dy/dt = -y³ com y(0) = 1, utilizando os métodos de Euler Explícito 
    e Euler Implícito. 
    
    Como o método implícito resulta em uma equação 
    algébrica não-linear para o passo futuro y_(i+1), utiliza-se o método 
    de Newton-Raphson para a convergência iterativa a cada passo de tempo.


    Analisa-se graficamente as aproximações numéricas 
    e a solução analítica exata y(t) = 1 / sqrt(2t + 1).
'''

import numpy as np 
import matplotlib.pyplot as plt 

def f(t, y):
    '''
    Definindo o pvi: dy/dt = -y^3; y(0) = 1
    '''
    v  = - (y ** 3)
    return v

def sol_analitica(t):
    '''
    Sol. analítica exata: y(t) = 1 / sqrt(2t + 1)
    '''
    y = 1.0 / np.sqrt(2 * t + 1)
    return y

def euler_explicito(f, t_final, dt, y0, k):
    '''
    Sol. numérica pelo Método de Euler Explícito para geração das curvas, em função dos passos que se dá, i.e. os deltihas (dt)
    '''
    t_array = np.arange(0, t_final + dt, dt) # N.B. --- dt garante o último ponto seja incluído, o np.array trunca no penúltimo!
    y_array = np.zeros(len(t_array))
    y_array[0] = y0
    
    for i in range(len(t_array) - 1):
        y_array[i + 1] = y_array[i] + dt * k * f(t_array[i], y_array[i])

    return t_array, y_array

t_e, y_e = euler_explicito(f, t_final=10e0, dt=0.25, y0=1e0, k=1e0)

def euler_implicito_nr(f, t_final, dt, y0, k):
    t_array = np.arange(0, t_final + dt, dt)
    y_array = np.zeros(len(t_array))
    y_array[0] = y0

    for i in range(len(t_array) - 1):
        erro_nr = 1e0
        y_new = y_array[i]

        while (erro_nr >= 1e-6):
            # f(t, y_new) já é negativo (-y^3), então para subtrair dt*f fica positivo:
            g = y_new - dt * k * f(t_array[i + 1], y_new) - y_array[i]
            
            # Derivada: g(y_new) = y_new + dt * k * y_new^3 - y_old
            dg = 1 + (3 * dt * k * (y_new ** 2)) 

            y_newer = y_new - g / dg
            erro_nr = abs(y_newer - y_new)
            y_new = y_newer

        # CORREÇÃO 1: Esta linha deve estar DENTRO do laço 'for' para salvar cada passo temporal
        y_array[i + 1] = y_new

    return t_array, y_array

# Mudança de t_final para 5e0 para alinhar com o xlim do gráfico
t_i, y_i = euler_implicito_nr(f, t_final=5e0, dt=0.25, y0=1e0, k=1e0)

ls_dt = [0.25, 0.55, 0.85, 1.05]
ls_t = np.arange(0, 5.1, 0.1) # Adicionado 5.1 para incluir o valor 5 no gráfico analítico
ls_sol = []

# CORREÇÃO 2: Alterado de 't' para 't_ponto' para não sobrescrever o vetor 't' do Euler Implícito
for t_ponto in ls_t:
    sol = sol_analitica(t_ponto)
    ls_sol.append(sol)


fig, axs = plt.subplots(2, 2, figsize=(10, 7)) # O figsize ajuda a não amassar os gráficos

for dt, ax in zip(ls_dt, axs.flatten()):

    t_e, y_e = euler_explicito(f, t_final=5e0, dt=dt, y0=1e0, k=1e0)
    t_i, y_i = euler_implicito_nr(f, t_final=5e0, dt=dt, y0=1e0, k=1e0)

    ax.plot(t_e, y_e, ls="--", marker="^", label="Euler Explícito")
    ax.plot(t_i, y_i, ls="--", marker="o", label="Euler Implícito N-R") 
    ax.plot(ls_t, ls_sol, color="black", label="Solução Analítica")
    ax.set_ylabel("y")
    ax.set_xlabel("t") 
    ax.set_title(f"h = {dt}")
    ax.set_xlim(0, 5)
    ax.legend() 

plt.tight_layout()
plt.savefig("Graphics/grafico_exercicio13.png", dpi=300, bbox_inches="tight")
plt.show()