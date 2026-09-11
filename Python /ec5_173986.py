"""
Disciplina: EQ 220 - Métodos Numéricos
Atividade: Exercício Computacional 5
Nome: Daniel Mussato Campiotti
Data: 09/09/2026
Descrição: Eliminação de Gauss e Substituição Reversa.
"""

# Matriz A (coeficientes) e Vetor b (termos independentes)

matriz = [
    [2.0, 1.0, -1.0],
    [-3.0, -1.0, 2.0],
    [-2.0, 1.0, 2.0]
]

vetor = [8.0, -11.0, -3.0] # i.e. "sol. do sistema"

n = len(vetor)

# eliminação de Gauss
h = 0
k = 0

while (h < n and k < n):
    i = h + 1
    while (i < n):
        f = matriz[i][k] / matriz[h][k]

        matriz[i][k] = 0.0  

        j = k + 1

        while (j < n):
            matriz[i][j] = matriz[i][j] - matriz[h][j] * f
            j += 1
            
        vetor[i] -= vetor[h] * f
        i += 1
        
    h += 1
    k += 1

# --- 2. PROCEDIMENTO SUBSTITUIÇÃO REVERSA ---
x = [0.0] * n
x[n - 1] = vetor[n - 1] / matriz[n - 1][n - 1]

i = n - 2
while (i >= 0):
    s = 0.0
    j = i + 1
    while (j < n):
        s = s + matriz[i][j] * x[j]
        j += 1
        
    x[i] = (vetor[i] - s) / matriz[i][i]
    i -= 1

# --- RESULTADOS ---
print("\nMatriz Triangularizada:")
for linha in matriz:
    print([round(v, 4) for v in linha])

print("\nVetor b Modificado:", [round(v, 4) for v in vetor])
print("\nSolução do Sistema (x):", [round(v, 4) for v in x])
print()