"""
Disciplina: EQ 220 - Métodos Numéricos
Atividade: Projeto Computacional 1
Nome:
    Daniel Mussato Campiotti
Data: 28/08/2026
Descrição: Interpolação polinomial para a modelagem do calor específico (à pressão constante) Cp do 
do sólido de Ni. 

    T (K) = 1500, 1900, 2300, 2700

    Cp (j.mol⁻¹.K⁻¹) = 29.4756, 32.160, 35.6409, 42.2840

Spline quadrático aos dados do ec9 utilizando SciPy.
    
"""

import os 
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline  

T = np.array([1500, 1900, 2300, 2700], dtype=float) # values of T (K)

Y = np.array([29.4756, 32.160, 35.6409, 42.2840], dtype=float) # values of Cp (j.mol⁻¹.K⁻¹)

n = len(T) # vector's lenght

spline_funcao = make_interp_spline(T, Y, k=2)  # k define o tipo de spline e.g. k=3 seria a Cúbica


ls_T = np.linspace(T[0], T[-1], 250)

# Avalia o spline diretamente passando o vetor de pontos desejados
ls_Cp = spline_funcao(ls_T) 


# Construção gráfico
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(ls_T, ls_Cp, ls="-", color="darkred", label="Spline Quadrática (SciPy)")
ax.scatter(T, Y, color="black", marker="o", s=60, zorder=5, label="Dados Experimentais (Ni)")

ax.set_ylabel("Cp ($J \cdot mol^{-1} \cdot K^{-1}$)")
ax.set_xlabel("T (K)")
ax.set_title("Modelagem do Calor Específico do Níquel Sólido")
ax.grid(True, linestyle=":", alpha=0.6)
ax.legend()

os.makedirs("Graphics", exist_ok=True)
plt.tight_layout()
plt.savefig("Graphics/projeto1_spline_scipy.png", dpi=300, bbox_inches="tight")
plt.show()
