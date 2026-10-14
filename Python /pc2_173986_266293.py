'''
Disciplina: Métodos Numéricos - EQ220
 Atividade: Projeto Computacional 2

 Nomes:
     - Daniel Mussato Campiotti (173986)
     - Vitoria Walderez Paulini (266293)

 Data: 12/10/2026

 Descrição: O trabalho consiste na resolução numérica de um sistema não linear
 de equações que descreve a transferência de calor em uma rede térmica.

 Para isso, implementou-se e comparou-se os métodos de Newton Multidimensional
 e de Broyden, utilizando como critério de convergência a norma do vetor de
 resíduos do sistema.

 O modelo considera mecanismos de transferência de calor por condução,
 convecção e radiação térmica, sendo a contribuição radiativa responsável
 pela não linearidade do sistema.

 Finalmente, analisou-se a convergência dos métodos, as soluções obtidas e a
 evolução dos resíduos ao longo das iterações.
'''

import numpy as np
import matplotlib.pyplot as plt


# Parâmetros físicos iniciais e constantes:

# Áreas[m²]
A = 1.0 
A_menor = A / 2

# Constantes
e = 0.8 # emissividade 
sigma = 5.670374419e-8 # constante de Stefan-Boltzmann 

# Parâmetros de troca de calor convectiva em [W m⁻² K⁻¹]
h_hot = 50.0 
h_cold = 25.0

# Temperaturas em [K] 
T_hot = 400.0 
T_cold = 300.0 




# Resistência termica
R1 = 0.02 / (10*A)
R2_1 = 0.015 / (15*A_menor) # leia-se R'2
R2_2 = 0.015 / (8*A_menor)  # leia-se R''2
R3 = 0.02 / (12*A)
R4_1 = 0.01 / (20*A_menor)
R4_2 = 0.01 / (10*A_menor)
R5 = 0.02 / (5*A)


# Definindo as funções necessaria para implementação do codigo.
def f(T):
    T0, T1, T2, T3, T4, T5 = T

    # Existem duas funções relacionadas às outras, i.e.:
    y = (T1 - T2) * (1/R2_1 + 1/R2_2)
    x = (T3 - T4) * (1/R4_1 + 1/R4_2)

    # Isolando-as, vem que:
    F0 = h_hot*A*(T_hot-T0) + e*sigma*A*(T_hot**4-T0**4) - (T0-T1)/R1
    F1 = (T0-T1)/R1 - y
    F2 = (T2-T3)/R3 - y
    F3 = (T2-T3)/R3 - x
    F4 = (T4-T5)/R5 - x
    F5 = (T4-T5)/R5 - h_cold*A*(T5-T_cold)

    # Finalmente, nosso vetor f é:
    return np.array([F0, F1, F2, F3, F4, F5])


# Aplicação do metodo de Newton Multidimensional, com o jacobiano inverso
def Jacobiano(T):
    T0, T1, T2, T3, T4, T5 = T

    # Derivadas parciais de cada função em cada uma das linhas
    y = (1/R2_1) + (1/R2_2)
    x = (1/R4_1) + (1/R4_2) # Constante

    J = np.array([
        [-h_hot*A - 4*e*sigma*A*T0**3 - 1/R1, 1/R1, 0.0, 0.0, 0.0, 0.0], # F0 em relação as temp
        [1/R1, -1/R1-y, y, 0.0, 0.0, 0.0],
        [0.0, -y, 1/R3+y, -1/R3, 0.0, 0.0],
        [0.0, 0.0, 1/R3, -1/R3-x, x, 0.0],
        [0.0, 0.0, 0.0, -x, 1/R5+x, -1/R5],
        [0.0, 0.0, 0.0, 0.0, 1/R5, -1/R5-h_cold*A]
    ])

    # Veja-se que se tem sobretudo uma matriz predominantemente diagonal, o que permitiria uma boa solução, também, por métodos iterativos!
    return J



def Inverso_J(T):
    '''
    Calcula-se o inverso do Jacobiano, para aplicar o algebrismo matricila da solução.
    '''
    T0, T1, T2, T3, T4, T5 = T

    jacobiano_inv = np.linalg.inv(Jacobiano(T))

    return jacobiano_inv


# Chute inicial
T0 = np.array([390.0, 380.0, 360.0, 340.0, 320.0, 310.0])
# É de se esperar que a temperatura decaia, por razões termodinâmicas! Daí a razão de escolher esse intervalo descrescente.

tolerancia = 1e-6


# Método de Newton Multidimensional
def Newton(f, T0, tol=1e-6):
    T = T0.copy()
    erro_newton = [np.linalg.norm(f(T), 2)]

    while (erro_newton[-1] >= tol):
        # Resolve o sistema linear J p = -f(T)
        # sem calcular explicitamente a inversa do Jacobiano.
        p = np.linalg.solve(Jacobiano(T), -f(T))
        T_novo = T + p
        erro_newton.append(np.linalg.norm(f(T_novo), 2))
        T = T_novo

    return T, erro_newton


# Aplicação do método de Newton
T, erro_newton = Newton(f, T0, tol=tolerancia)
print("Solução pelo método de Newton:")
print(T)


def Jacobiano_perturbado(f, T0, epsilon=1e-6):
    '''
    Aproximação numérica do Jacobiano por uma perturbação epsilon.

    N.B. --- Percebeu-se que, o tentar esboçar alguma matriz como np.eye(n) * T_hot, prejudicou sobremaneira a convergência do método de Broyden.
    A fim, portanto, de evitar isso, causar-se-á através dessa perturbação infinetesimal uma composição de um Jacobiano grosseiro (H), cujo uso será
    destinado à Broyden. Broyden exige que se tenha um chute inicial de uma matriz Jacobiana.

    Não sabemos ao certo se essa confusão de métodos é, contudo, coerente. Entre o uso de expressões bizarras para o chute inicial do Jacobiano este
    foi a único chute que queremos crer coerente; tanto ao fenômeno, quanto ao método numérico. Escolher de partida o Jacobiano analítico seria um erro
    metodológico grave e escolher uma matriz qq. preenchida por parâmetros incertos e conduzidos pelos valores vistos nesse Jacobiano anlítico seria ainda falso.
    '''
    n = len(T0)
    H = np.zeros((n, n)) # Jacobiano perturbado (H)
    f_inic = f(T0)

    for k in range(n):
        T_pert = T0.copy()

        # perturbação infinetsimal por epsilon na Temp. j
        T_pert[k] += epsilon   # perturbação na temperatura
        f_pert = f(T_pert) # célculo após a perturbação

        # alocação por coluna no jacobiano perturbado
        H[:, k] = (f_pert - f_inic)/epsilon

    return H

# Método de Broyden para resolver o sistema de equações não lineares
def Broyden(f, T0, tol=1e-6):
    '''
    Resolve um sistema não linear pelo método de Broyden.

    Referência:
        BROYDEN, C. G. (1965). A class of methods for solving
        nonlinear simultaneous equations.
        https://doi.org/10.1090/S0025-5718-1965-0198670-6
    '''
    n = len(T0)
    T = T0.copy()

    # Aproximação inicial de B
    B = Jacobiano_perturbado(f, T0, epsilon=1e-6)

    erros = []

    # Resíduo inicial do sistema
    f_k = f(T)
    erro = np.linalg.norm(f_k)
    erros.append(erro)

    while (erro >= tol):
        # Resolve B p_k = -f(T)
        p_k = np.linalg.solve(B, -f_k)
        T_novo = T + p_k
        f_novo = f(T_novo)

        # Variação do vetor de resíduos
        y_k = f_novo - f_k

        # Produto da matriz aproximada pelo passo
        Bp = B @ p_k
        denominador = np.dot(p_k, p_k)

        if denominador == 0:
            raise ValueError("Passo nulo.")

        # Atualização de Broyden o
        for i in range(n):
            for j in range(n):
                B[i, j] += (y_k[i] - Bp[i]) * p_k[j] / denominador

        # Erro como norma do vetor de resíduos
        erro = np.linalg.norm(f_novo)
        erros.append(erro)

        # Atualização das temperaturas e dos resíduos
        T = T_novo
        f_k = f_novo

    return T, B, erros


# Aplicação do método de Broyden
T_broyden, J_broyden, erro_broyden = Broyden(f, T0, tol=tolerancia)


# Comparação dos resultados
# Jacobiano analítico avaliado na solução de Broyden
J_analitico = Jacobiano(T_broyden)

print("\nSolução pelo método de Newton:")
print(T)
print("\nResíduo final de Newton:")
print(np.linalg.norm(f(T), 2))
print("\nNúmero de iterações de Newton:")
print(len(erro_newton) - 1)

print("\nSolução pelo método de Broyden:")
print(T_broyden)
print("\nResíduo final de Broyden:")
print(np.linalg.norm(f(T_broyden), 2))
print("\nNúmero de iterações de Broyden:")
print(len(erro_broyden) - 1)

print("\nJacobiano analítico na solução de Broyden:")
print(J_analitico)
print("\nJacobiano aproximado pelo método de Broyden:")
print(J_broyden)
print("\nNorma da diferença entre os jacobianos:")
print(np.linalg.norm(J_analitico - J_broyden))
print("\nNorma da diferença entre as soluções:")
print(np.linalg.norm(T - T_broyden))


# Gráfico do erro em função do número de iterações
# Índices das iterações, incluindo a estimativa inicial (iteração 0)
iteracoes_newton = np.arange(len(erro_newton))
iteracoes_broyden = np.arange(len(erro_broyden))

plt.figure(figsize=(9, 6))
plt.semilogy(iteracoes_newton, erro_newton, marker='o', markersize=4,
             label='Newton', color='purple')
plt.semilogy(iteracoes_broyden, erro_broyden, marker='s', markersize=4,
             label='Broyden', color='blue')
plt.axhline(tolerancia, linestyle='--', label=f'Tolerância = {tolerancia:.0e}')
plt.xlabel('Número de iterações')
plt.ylabel(r'Erro: $\|f(T_k)\|$ (W)')
plt.title('Convergência dos métodos de Newton e Broyden')
plt.grid(True, which='both', linestyle=':')
plt.legend()
plt.tight_layout()
plt.show()


# Análise de linearidade tolerância e núemro de iteração
ls_tolerancia = np.logspace(0, -10, num=10)
ls_erro_newton = []
ls_erro_broyden = []

for tolerancia_atual in ls_tolerancia:
    T_n, erro_n = Newton(f, T0, tol=tolerancia_atual)
    T_b, B_b, erro_b = Broyden(f, T0, tol=tolerancia_atual)
    ls_erro_newton.append(len(erro_n) - 1)
    ls_erro_broyden.append(len(erro_b) - 1)

print("Iterações de Newton:", ls_erro_newton)
print("Iterações de Broyden:", ls_erro_broyden)


# print gráfico N_iterações vs. tolerancia
plt.figure(figsize=(8, 5))
plt.semilogx(ls_tolerancia, ls_erro_newton, 'o-',
             label='Newton-Raphson', color='blue')
plt.semilogx(ls_tolerancia, ls_erro_broyden, 's-',
             label='Broyden', color='purple')
plt.xlabel('Tolerância')
plt.ylabel('Número de iterações')
plt.title('Influência da tolerância no número de iterações')
plt.grid(True)
plt.legend()
plt.gca().invert_xaxis()
plt.tight_layout()
plt.show()

'''
Os resultados obtidos mostram que os métodos de Newton multidimensional e de Broyden convergiram para praticamente a mesma solução do sistema não linear, apresentando resíduos finais da ordem de (10⁻¹¹). 
Para o chute inicial e a tolerância considerados, o método de Newton necessitou de três iterações, enquanto o método de Broyden necessitou de quatro.

Essa diferença é compatível com as características dos métodos: Newton utiliza o Jacobiano analítico diretamente atualizado a cada iteração, 
enquanto Broyden atualiza uma aproximação do Jacobiano com base nas variações do vetor de resíduos.
Embora Broyden possa reduzir o custo associado ao recálculo do Jacobiano em problemas maiores, 
seu desempenho depende, sobremaneira, da aproximação inicial utilizada e do comportamento do sistema.

Na análise da influência da tolerância, observou-se que o número de iterações não varia linearmente, apresentando um comportamento em degraus.
Isso ocorre porque o número de iterações é inteiro e depende de quantas etapas são necessárias para que o resíduo fique abaixo de cada tolerância.
Os resultados, portanto, indicam que ambos os métodos foram eficazes para a rede térmica estudada, com vantagem em número de iterações para Newton neste caso específico.
'''
