"""
Disciplina: EQ 220 - Métodos Numéricos
Atividade: Exercício Computacional 11
Nome:
    Daniel Mussato Campiotti
Data: 28/08/2026
Descrição:
    Ajuste linear da cinética de crescimento bacteriano comparando o
    MMQ analítico clássico com a minimização numérica via SciPy.
"""

import numpy as np
from scipy.optimize import minimize

t = np.linspace(3, 11, num=9)
X = np.log(np.array([0.03, 0.06, 0.16, 0.43, 0.35, 1.13, 2.58, 2.71, 5.34]))
n = len(t)

def mmq(t, X):
    '''Cálculo direto do mmq.'''
    a1 = (np.mean(t * X) - np.mean(X) * np.mean(t)) / (np.mean(t ** 2) - (np.mean(t)) ** 2)
    a0 = np.mean(X) - a1 * np.mean(t)
    return a0, a1

a0, a1 = mmq(t, X)

def R2(t, X, a0, a1):
    '''Cálculo do R² (Ajustado para consistência estatística)'''
    # Soma dos Quadrados dos Resíduos (Sem a divisão por n para bater com SStot)
    ss_res = np.sum((X - a0 - a1 * t) ** 2)

    # Soma dos Quadrados Totais
    ss_tot = np.sum((X - np.mean(X)) ** 2)
    
    r2 = 1 - (ss_res / ss_tot)
    return r2

r2 = R2(t, X, a0, a1)

print("--- MÉTODO ANALÍTICO ---")
print(f"a₀ = {a0:.6f}")
print(f"a₁ = {a1:.6f}")
print(f"Equação: ln[X] = {a0:.4f} + {a1:.4f} ⋅ t")
print(f"R²:      {r2:.6f}")

# Método de ajuste por scipy.optimize.minimize

def custo(coef, t, X):
    a0_min = coef[0]
    a1_min = coef[1]
    C = np.sum((X - a0_min - a1_min * t)**2)
    return C

chute = [0.0, 0.0]

# min(C(a0,a1)) é definido como:
res_otim = minimize(custo, chute, args=(t, X))

a0_min, a1_min = res_otim.x
r2_min = R2(t, X, a0_min, a1_min)

print("\n--- MÉTODO NUMÉRICO (SCIPY) ---")
print(f"a₀ = {a0_min:.6f}")
print(f"a₁ = {a1_min:.6f}")
print(f"Equação: ln[X] = {a0_min:.4f} + {a1_min:.4f} ⋅ t")
print(f"R²:      {r2_min:.6f}")

# Cálculo das diferenças absolutas (Δ) para comparação da convergência
print("\n--- ANÁLISE DE DIFERENÇAS (Δ) ---")
print(f"|Δa₀|: {abs(a0 - a0_min):.2e}")
print(f"|Δa₁|: {abs(a1 - a1_min):.2e}")
print(f"|ΔR²|: {abs(r2 - r2_min):.2e}")
