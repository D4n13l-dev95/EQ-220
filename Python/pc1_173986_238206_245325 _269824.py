"""
Disciplina: EQ 220 - Métodos Numéricos

Atividade: Projeto Computacional 1

Nome:

    Daniel Mussato Campiotti (173986)

    Jeferson dos Santos Teixeira (238206)

    Miguel Teodoro Souza (269824)

    Pedro Fernandes dos Santos (245325)

Data: 21/08/2026

Descrição: Determinação da pressão de saturação da água através da
    equação de van der Waals. As raízes da equação cúbica em Z são determinadas 
    pelo método de Newton-Raphson. 
    A escolha do P_sat deu-se pela aplicação do regula falsi
    a condição de equilíbrio líquido-vapor é:
        phi_liq = phi_vap
    O cotejo com os dados experimentais do NIST, foi através do plot dum gráfico.
"""

import numpy as np

import matplotlib.pyplot as plt

# funções da equação de van der Waals

def F(Z, A, B):

    Y = Z**3 - (1 + B) * Z**2 + A * Z - A * B

    return Y

def dF(Z, A, B):

    dY = 3 * Z**2 - 2 * (1 + B) * Z + A

    return dY

def ln_phi(Z, A, B):

    # calcula o logaritmo do coeficiente de fugacidade

    phi = Z - 1 - np.log(Z - B) - A / Z

    return phi

# constantes

R = 8.314462       # cm3 MPa K^-1 mol^-1

Tc = 647.1         # temperatura crítica [K]

Pc = 22.06         # pressão crítica [MPa]

# parâmetros a e b da equação de van der Waals

a = 27 * R**2 * Tc**2 / (64 * Pc)

b = R * Tc / (8 * Pc)

# dados experimentais do NIST

# temperatura [°C] vs. pressão [kPa]

nist_data = {
    20: 2.339,

    25: 3.169,

    30: 4.247,

    35: 5.629,

    40: 7.385,

    45: 9.595,

    50: 12.35,

    55: 15.76,

    60: 19.95,

    65: 25.04,

    70: 31.20,

    75: 38.58,

    80: 47.41,

    85: 57.87,

    90: 70.18,

    95: 84.61,

    100: 101.42
}

# tolerâncias utilizadas nos métodos numéricos

tolerancia_Z = 1e-12

tolerancia_P = 1e-10

max_iteracoes = 200

# temperaturas utilizadas no cálculo

lista_T = np.arange(20, 105, 5)

# listas para armazenar os resultados

lista_Psat_vdw = []

lista_Psat_nist = []

# cálculo para cada temperatura

for T_celsius in lista_T:

    # conversão da temperatura de °C para K

    T = T_celsius + 273.15

    # procura da pressão de saturação
    # usamos vários valores de pressão para encontrar uma mudança
    # de sinal na diferença entre as fugacidades

    lista_P = np.geomspace(1e-5, 5.0, 10000)

    P_anterior = None

    g_anterior = None

    for P in lista_P:

        # cálculo dos parâmetros adimensionais A e B

        A = a * P / (R**2 * T**2)

        B = b * P / (R * T)

        # raiz líquida por Newton-Raphson
        
        # chute inicial ligeiramente maior que B,
        # pois necessário ter Z > B para calcular ln(Z - B)

        Z = B + 1e-4 # mínimo Z para não dar problema no ln

        contador = 0

        erro = 1.0

        while (erro >= tolerancia_Z and contador < max_iteracoes):

            derivada = dF(Z, A, B)

            # evita divisão por uma derivada muito próxima de zero

            if (abs(derivada) < 1e-12):

                Z = None

                break

            # fórmula do método de Newton-Raphson

            Z_novo = Z - F(Z, A, B) / derivada

            # verifica-se se o novo valor de Z continua válido

            if (Z_novo <= B or not np.isfinite(Z_novo)):

                Z = None

                break

            # erro utilizado para verificar a convergência

            erro = abs(Z_novo - Z)

            Z = Z_novo

            contador += 1

        Z_liq = Z
        
        # raiz de vapor por Newton-Raphson
        
        # para a fase vapor, Z costuma estar próximo de 1

        Z = 1.0 # máximo valor de Z

        contador = 0

        erro = 1.0

        while (erro >= tolerancia_Z and contador < max_iteracoes):

            derivada = dF(Z, A, B)

            # evita divisão por uma derivada muito próxima de zero

            if (abs(derivada) < 1e-12):

                Z = None
                break

            # fórmula do método de Newton-Raphson

            Z_novo = Z - F(Z, A, B) / derivada

            # verifica se o novo valor de Z continua válido

            if (Z_novo <= B or not np.isfinite(Z_novo)):
                Z = None
                break

            erro = abs(Z_novo - Z)

            Z = Z_novo

            contador += 1

        Z_vap = Z

        # verifica se as duas raízes foram encontradas

        if Z_liq == None or Z_vap == None:

            continue

        # garante que Z_liq seja a menor das duas raízes encontradas

        if (Z_liq > Z_vap):

            Z_liq, Z_vap = Z_vap, Z_liq

        # cálculo das fugacidades
        
        ln_phi_liq = ln_phi(Z_liq, A, B)

        ln_phi_vap = ln_phi(Z_vap, A, B)

        # diferença entre as fugacidades
        # no eq. líquido-vapor, essa diferença deve ser zero

        g_atual = ln_phi_liq - ln_phi_vap
        
        # procura de uma mudança de sinal
        
        # quando g muda de sinal, sabemos que a pressão de saturação
        # está entre a pressão anterior e a pressão atual

        if (P_anterior is not None):

            if (g_anterior * g_atual < 0):

                # encontramos o intervalo que contém a solução

                P_inferior = P_anterior

                P_superior = P

                break

        P_anterior = P

        g_anterior = g_atual

    # método da regula falsi para P_sat

    erro = 1.0

    contador = 0

    P_anterior = P_superior

    while (erro >= tolerancia_P and contador < max_iteracoes):

        # cálculo de A e B no limite inferior

        A_inferior = a * P_inferior/(R**2*T**2)

        B_inferior = b * P_inferior/(R*T)

        # raiz líquida no limite inferior

        Z = B_inferior + 1e-4

        erro_Z = 1.0

        contador_Z = 0

        while erro_Z >= tolerancia_Z and contador_Z < max_iteracoes:

            derivada = dF(Z, A_inferior, B_inferior)

            Z_novo = Z - F(Z, A_inferior, B_inferior) / derivada

            erro_Z = abs(Z_novo - Z)

            Z = Z_novo

            contador_Z += 1

        Z_liq_inferior = Z

        # raiz de vapor no limite inferior

        Z = 1.0

        erro_Z = 1.0

        contador_Z = 0

        while erro_Z >= tolerancia_Z and contador_Z < max_iteracoes:

            derivada = dF(Z, A_inferior, B_inferior)

            Z_novo = Z - F(Z, A_inferior, B_inferior) / derivada

            erro_Z = abs(Z_novo - Z)

            Z = Z_novo

            contador_Z += 1

        Z_vap_inferior = Z

        # diferença de fugacidade no limite inferior

        g_inferior = ( 
            ln_phi(Z_liq_inferior, A_inferior, B_inferior)
            -
            ln_phi(Z_vap_inferior, A_inferior, B_inferior)
            )

        # cálculo de A e B no limite superior

        A_superior = a * P_superior/(R**2*T**2)

        B_superior = b * P_superior/(R*T)

        # raiz líquida no limite superior

        Z = B_superior + 1e-4

        erro_Z = 1.0

        contador_Z = 0

        while erro_Z >= tolerancia_Z and contador_Z < max_iteracoes:

            derivada = dF(Z, A_superior, B_superior)

            Z_novo = Z - F(Z, A_superior, B_superior) / derivada

            erro_Z = abs(Z_novo - Z)

            Z = Z_novo

            contador_Z += 1

        Z_liq_superior = Z

        # raiz de vapor no limite superior

        Z = 1.0

        erro_Z = 1.0

        contador_Z = 0

        while erro_Z >= tolerancia_Z and contador_Z < max_iteracoes:

            derivada = dF(Z, A_superior, B_superior)

            Z_novo = Z - F(Z, A_superior, B_superior) / derivada

            erro_Z = abs(Z_novo - Z)

            Z = Z_novo

            contador_Z += 1

        Z_vap_superior = Z

        # diferença de fugacidade no limite superior

        g_superior = (
            ln_phi(Z_liq_superior, A_superior, B_superior)
            -
            ln_phi(Z_vap_superior, A_superior, B_superior)
        )

        # calcula o novo valor de pressão pela regula falsi

        P_novo = ( P_inferior * g_superior - P_superior * g_inferior) / (g_superior - g_inferior)

        # cálculo de A e B no novo valor de pressão

        A = a * P_novo / (R**2 * T**2)

        B = b * P_novo / (R * T)

        # raiz líquida

        Z = B + 1e-4

        erro_Z = 1.0

        contador_Z = 0

        while erro_Z >= tolerancia_Z and contador_Z < max_iteracoes:

            derivada = dF(Z, A, B)

            Z_novo = Z - F(Z, A, B) / derivada

            erro_Z = abs(Z_novo - Z)

            Z = Z_novo

            contador_Z += 1

        Z_liq = Z

        # raiz de vapor

        Z = 1.0

        erro_Z = 1.0

        contador_Z = 0

        while erro_Z >= tolerancia_Z and contador_Z < max_iteracoes:

            derivada = dF(Z, A, B)

            Z_novo = Z - F(Z, A, B) / derivada

            erro_Z = abs(Z_novo - Z)

            Z = Z_novo

            contador_Z += 1

        Z_vap = Z

        # diferença de fugacidade no novo valor de pressão

        g_novo = (
            ln_phi(Z_liq, A, B)
            -
            ln_phi(Z_vap, A, B)
        )

        # verifica em qual lado da raiz está o novo ponto

        if g_inferior * g_novo < 0:

            P_superior = P_novo

        else:

            P_inferior = P_novo

        # erro entre duas aproximações consecutivas

        erro = abs(P_novo - P_anterior)

        P_anterior = P_novo

        contador += 1

    # pressão de saturação calculada pela equação de van der Waals

    P_sat = P_novo

    # conversão de MPa para kPa

    P_sat_kPa = P_sat * 1000

    # valor experimental correspondente à temperatura

    P_nist = nist_data[T_celsius]

    # cálculo do erro percentual em relação ao valor experimental

    erro_percentual = (abs(P_sat_kPa - P_nist)/P_nist) * 100

    # armazena os resultados para o gráfico e para o R²

    lista_Psat_vdw.append(P_sat_kPa)

    lista_Psat_nist.append(P_nist)

    # mostra os resultados para cada temperatura

    print(
        f"T = {T_celsius:3.0f} °C | "

        f"P_vdw = {P_sat_kPa:8.3f} kPa | "

        f"P_NIST = {P_nist:8.3f} kPa | "

        f"erro = {erro_percentual:7.2f} %"
    )

# Interlúdio e Discussão dos resultados

# Os resultados mostram que a pressão de saturação aumenta com a
# temperatura, acompanhando a tendência observada nos dados do NIST.
# entretanto, os valores calculados pela equação de van der Waals
# apresentam desvios consideráveis em relação aos valores experimentais.
# Acreditamos, sobretudo, que essa diferença releaciona-se às simplificações da equação de
# van der Waals para representar as interações intermoleculares da água; i.e., o modelo matemático utilizado (polinomial c/ Z adimensional)
# não é tão próximo da realidade quanto se gostraria. Afinal, sendo a água o objeto de nosso
# modelo, sempre há desvios devido à natureza complexa desta.

# Apesar dos desvios quantitativos, o método reproduz o comportamento
# crescente da pressão de saturação com a temperatura.
# Do ponto de vista numérico, o método de Newton-Raphson foi utilizado
# para determinar as raízes da equação cúbica em Z, enquanto a regula
# falsi foi utilizada para determinar a pressão de saturação a partir
# da condição de igualdade das fugacidades entre as fases.

# N.B. --- A equação cúbica pode apresentar três raízes reais. No código, utilizamos dois chutes iniciais 
# distintos no método de Newton-Raphson:  um próximo de B, buscando a menor raiz associada à fase líquida,
# e outro próximo de Z = 1, buscando a maior raiz associada à fase vapor. 
# Assim, a raiz intermediária não é utilizada no cálculo das fugacidades. Quer dizer, o que quisemos 
# foi, em essência tentar evitá-la. A  geometria da eq de 3° grau que aqui se mostra, perimite esta razoável ideia,
# que queremos crer ter sido suficiente para Z_liq e Z_vap.

# gráfico

plt.figure(figsize=(8, 6))

plt.plot(
    lista_T,
    lista_Psat_vdw,
    'o-',
    label='van der Waals',
    color="purple"
)

plt.plot(
    lista_T,
    lista_Psat_nist,
    's--',
    label='NIST',
    color="blue"
)

plt.xlabel('Temperatura (°C)')
plt.ylabel('Pressão de saturação (kPa)')
plt.title('Pressão de saturação da água')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
