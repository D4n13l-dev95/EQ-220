'''
Disciplina: EQ 220 - Métodos Numéricos
Atividade: Exercício Computacional 12
Nome:
    Daniel Mussato Campiotti
Data: 27/10/2026
Descrição: Resolução do seguinte problema da valor inicial:
    dy / dt = -y³
em que,
    y(0) = 1.
    Pelo Método de Euler Explícito.
'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

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


def euler_curvas(f, t_final, dt, y0):
    '''
    Sol. numérica pelo Método de Euler Explícito para geração das curvas, em função dos passos que se dá, i.e. os deltihas (dt)
    '''

    t_array = np.arange(0, t_final + dt, dt) # N.B. --- dt garante o último ponto seja incluído, o np.array trunca no penúltimo!

    y_array = np.zeros(len(t_array))

    y_array[0] = y0
    
    for i in range(len(t_array) - 1):

        y_array[i + 1] = y_array[i] + dt * f(t_array[i], y_array[i])

    return t_array, y_array

def R_2(y_real, y_pred):
    '''
    Cálculo do R² comparando pontos nos mesmos instantes de tempo
    '''

    y_media = np.mean(y_real)

    ss_tot = np.sum((y_real - y_media) ** 2)

    ss_res = np.sum((y_real - y_pred) ** 2)

    r_2 = 1 - (ss_res / ss_tot)

    return r_2

# Para a "simulação"
t_final = 10.0

# Execução do método de Euler para diferentes tamanhos de passo (dt)
t_grande, y_grande = euler_curvas(f, t_final, dt=1.05, y0=1.0)
t_medio, y_medio   = euler_curvas(f, t_final, dt=0.8, y0=1.0)
t_suave, y_suave   = euler_curvas(f, t_final, dt=0.1, y0=1.0)

# Malha contínua para plotagem da solução analítica
t_analitico = np.linspace(0, t_final, 1000)
y_analitico = sol_analitica(t_analitico)

# Cálculo correto do R² para cada passo
r2_grande = R_2(sol_analitica(t_grande), y_grande)
r2_medio  = R_2(sol_analitica(t_medio), y_medio)
r2_suave  = R_2(sol_analitica(t_suave), y_suave)

print(f"R² (dt = 1.05): {r2_grande:.6f}")
print(f"R² (dt = 0.80): {r2_medio:.6f}")
print(f"R² (dt = 0.10): {r2_suave:.6f}")

# --- Construção do Gráfico ---
fig, ax = plt.subplots(figsize=(10, 6), dpi=120)

ax.plot(t_analitico, y_analitico, label="Solução Analítica (Exata)", color="black", linewidth=2.5, zorder=3)
ax.plot(t_grande, y_grande, 'o--', label=f"Euler (dt = 1.05) | R² = {r2_grande:.4f}", color="purple", alpha=0.85)
ax.plot(t_medio, y_medio, 's--', label=f"Euler (dt = 0.80) | R² = {r2_medio:.4f}", color="blue", alpha=0.85)
ax.plot(t_suave, y_suave, '^:', label=f"Euler (dt = 0.10) | R² = {r2_suave:.4f}", color="green", alpha=0.85)

ax.set_title("Convergência do Método de Euler para $dy/dt = -y^3$", fontsize=13, fontweight='bold', pad=12)
ax.set_xlabel("Tempo ($t$)", fontsize=11)
ax.set_ylabel("$y(t)$", fontsize=11)
ax.set_xlim(-0.1, 10.1)
ax.set_ylim(-0.5, 1.05)
ax.grid(True, linestyle="--", alpha=0.5)
ax.legend(loc="upper right", frameon=True, fontsize=10)

plt.tight_layout()
plt.show()

# Relação entre dt e R² ?
dt_list = np.arange(0.05, 1.05, 0.05)

dados = []

for dt in dt_list:
   
    t, y_pred = euler_curvas(f, t_final=2, dt=dt)
    

    y_exato = sol_analitica(t)
    
    
    r2 = R_2(y_exato, y_pred)
    
    dados.append({"dt": round(dt, 3), "R2": r2})

#  Conversão para DataFrame Pandas e exibe a tabela
df = pd.DataFrame(dados)
print(df.to_string(index=False))

# Plotagem via Pandas
ax = df.plot(x="dt", y="R2", marker="o", figsize=(8, 4.5), grid=True)
ax.set_title("Variação do $R^2$ em Função de $dt$")
ax.set_xlabel("Passo de Integração ($dt$)")
ax.set_ylabel("Coeficiente $R^2$")
ax.set_ylim(-1,5)
plt.show()