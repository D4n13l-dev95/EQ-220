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

Primeiramente, é necessário calcular os polinômios de Lagrange para cada ponto de dados
Segundamente, é necessário plotar o polinômio de Lagrange e os pontos de dados para verificar a interpolação.
E, finalmente, coparar com o Método Direto de interpolação polinomial. 

"""

# listagem inicial dos pontos 

import numpy as np
import matplotlib.pyplot as plt

# initial data point and definitions used throughout the entire code

X = np.array([1500, 1900, 2300, 2700], dtype=float) # values of T (K)

Y = np.array([29.4756, 32.160, 35.6409, 42.2840], dtype=float) # values of Cp (j.mol⁻¹.K⁻¹)

n = len(X) # vector's lenght

b = [] # coefficientes of the polynomial interpolation (go to: line 107)



# iterations by def

def lagrange_interpolation(x_var, X, Y):
    """
    calculate the Lagrange interpolation polynomial.

    the Lagrange basis polynomials are calculated according to

        li = product[(x - X[j]) / (X[i] - X[j])]

    and the interpolation polynomial is calculated by

        L(x) = sum[Y[i] * li]
    """

    li = np.ones(n) # one's array for the Lagrange polynomials

    # calculate products for each Lagrange polynomial

    for i in range(n):

        for j in range(n):

            if (j != i):

                li[i] *= (x_var - X[j]) / (X[i] - X[j])

    # calculate the Lagrange interpolation polynomial
    # L = \sum_{i=0}^{n-1} Y[i] * li[i]

    L = np.dot(Y, li)

    return L

# N.B. - Notice that the polinomial expression is implicity. 
# So for whatever T necessary, one can use the definition above to calculate the value of Cp.
# Because li is a function of x_var, the value of Cp, so, can be calculated for any value of T, 
# not only for the points of data. 


# Direct Method of Polynomial Interpolation

# V . a = f, in which V is Vandermonde's matrix

V = np.vander(X, N=n) # Vandermonde's matrix (apparatley, Numpy already has its function)

def solve_direct_method(V, Y):
    """
    solve the system of equations V . a = Y for the coefficients a.
    """
    a = np.linalg.solve(V, Y)

    return a

# printing and calculating results

a = solve_direct_method(V, Y)

print("\nVetor solução (Método Direto):", a)

for i in range(n):

    print(f"a_direct_method[{i}] = {a[i]:.10f}")



# comparison with the Lagrange interpolation polynomial.
# discovering the coefficients of each Lagrange basis polynomial

def li(X):
    """
    calculate the coefficients of the Lagrange basis polynomials.
    """

    b = np.zeros((n, n))

    for i in range(n):

        # polynomial starts as 1
        coeff = np.array([1.0])

        for j in range(n):

            if (j != i):

                # multiply the polynomial by (x - X[j])
                new_coeff = np.zeros(len(coeff) + 1)

                for k in range(len(coeff)):
                    new_coeff[k] += coeff[k]
                    new_coeff[k + 1] -= coeff[k] * X[j]

                coeff = new_coeff

                # denominator of the Lagrange basis polynomial
                coeff /= (X[i] - X[j])

        b[i, :] = coeff

    return b

# printing and calcultaing coefficientes

b = li(X)

b_lagrange = np.dot(Y, b)

print("\nVetor solução (Lagrange):", b_lagrange)
for i in range(n):
    print(f"b_lagrange[{i}] = {b_lagrange[i]:.10f}")

# checking any distortion in data by norm of the difference of vectors

norm =np.linalg.norm(b_lagrange -a)
print(f"\nnorma(b_lagrange - a_direct_method) = {norm:.13f}")
print()




# graphical representation of the interpolation polynomial and the data points

x_plot = np.linspace(np.min(X), np.max(X), 100) # 100 points between the minimum and maximum of X

y_plot = np.array([
    lagrange_interpolation(x, X, Y)
    for x in x_plot
])


plt.scatter(X, 
            Y, 
            color='purple', 
            marker='o', 
            s=100, 
            label=r"Pontos Experimentais $(T, C_p)$")

plt.plot(x_plot,
         y_plot,
         color='black', 
         linestyle='-', 
         linewidth=2, 
         label=r"Polinômio de Lagrange $L(x) = \sum Y_i \cdot l_i(x)$")


plt.title(r"Interpolação Polinomial de Lagrange: $C_p$ do sólido de Ni", fontsize=12)

plt.xlabel(r'$T / K$', fontsize=11)

plt.ylabel(r'$C_p / (j \cdot mol^{-1} \cdot K^{-1})$', fontsize=11)

plt.legend(loc='upper left')

plt.grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()

plt.savefig("Graphics/interpolacao_lagrange.png", dpi=300)

plt.show()