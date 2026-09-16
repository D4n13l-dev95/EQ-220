"""
Disciplina: EQ 220 - Métodos Numéricos

Atividade: Projeto Computacional 1

Nome:

    Daniel Mussato Campiotti

    Jeferson

    Miguel Teodoro Souza

    Pedro

Data: 21/08/2026

Descrição:

    determinação da pressão de saturação da água através da

    equação de van der Waals.

    as raízes da equação cúbica em Z são determinadas pelo

    método de Newton-Raphson.

    a condição de equilíbrio líquido-vapor é:

        phi_liq = phi_vap

    os resultados são comparados com dados experimentais do NIST.

"""

import numpy as np

import matplotlib.pyplot as plt


# funções da equação de van der Waals

def F(Z, A, B):

    Y = Z**3 - (1 + B)*Z**2 + A*Z - A*B

    return Y


def dF(Z, A, B):

    dY = 3*Z**2 - 2*(1 + B)*Z + A

    return dY


def ln_phi(Z, A, B):

    # calcula o logaritmo do coeficiente de fugacidade

    phi = Z - 1 - np.log(Z - B) - A/Z

    return phi


# constantes

R = 8.314462       # cm3 MPa K^-1 mol^-1

Tc = 647.1         # temperatura crítica [K]

Pc = 22.06         # pressão crítica [MPa]


# parâmetros a e b da equação de van der Waals

a = 27*R**2*Tc**2/(64*Pc)

b = R*Tc/(8*Pc)


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

tolerancia_Z = 1e-10

tolerancia_P = 1e-8

max_iteracoes = 100


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

    lista_P = np.geomspace(1e-5, 5.0, 300)

    P_anterior = None

    g_anterior = None


    for P in lista_P:

        # cálculo dos parâmetros adimensionais A e B

        A = a*P/(R**2*T**2)

        B = b*P/(R*T)


        
        # raiz líquida por Newton-Raphson
        

        # chute inicial ligeiramente maior que B,
        # pois necessário ter Z > B para calcular ln(Z - B)

        Z = B + 1e-4

        contador = 0

        erro = 1.0


        while (erro >= tolerancia_Z and contador < max_iteracoes):

            derivada = dF(Z, A, B)

            # evita divisão por uma derivada muito próxima de zero

            if (abs(derivada) < 1e-12):

                Z = None

                break


            # fórmula do método de Newton-Raphson

            Z_novo = Z - F(Z, A, B)/derivada


            # verifica se o novo valor de Z continua válido

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

        Z = 1.0

        contador = 0

        erro = 1.0


        while (erro >= tolerancia_Z and contador < max_iteracoes):

            derivada = dF(Z, A, B)

            # evita divisão por uma derivada muito próxima de zero

            if (abs(derivada) < 1e-12):

                Z = None

                break


            # fórmula do método de Newton-Raphson

            Z_novo = Z - F(Z, A, B)/derivada


            # verifica se o novo valor de Z continua válido

            if (Z_novo <= B or not np.isfinite(Z_novo)):

                Z = None

                break


            erro = abs(Z_novo - Z)

            Z = Z_novo

            contador += 1


        Z_vap = Z


        # verifica se as duas raízes foram encontradas

        if Z_liq is None or Z_vap is None:

            continue


        # garante que Z_liq seja a menor das duas raízes encontradas

        if Z_liq > Z_vap:

            Z_liq, Z_vap = Z_vap, Z_liq


        # ---------------------------------
        # cálculo das fugacidades
        # ---------------------------------

        ln_phi_liq = ln_phi(Z_liq, A, B)

        ln_phi_vap = ln_phi(Z_vap, A, B)


        # diferença entre as fugacidades

        # no equilíbrio líquido-vapor, essa diferença deve ser zero

        g_atual = ln_phi_liq - ln_phi_vap


        
        # procura de uma mudança de sinal
        

        # quando g muda de sinal, sabemos que a pressão de saturação
        # está entre a pressão anterior e a pressão atual

        if (P_anterior is not None):

            if (g_anterior*g_atual < 0):

                # encontramos o intervalo que contém a solução

                P_inferior = P_anterior

                P_superior = P

                break


        P_anterior = P

        g_anterior = g_atual


    
    # método da bisseção para P_sat


    erro = 1.0

    contador = 0


    while (erro >= tolerancia_P and contador < max_iteracoes):

        # calcula o ponto médio do intervalo de pressão

        P_meio = (P_inferior + P_superior)/2


        # cálculo de A e B no ponto médio

        A = a*P_meio/(R**2*T**2)

        B = b*P_meio/(R*T)


    
        # raiz líquida
    

        Z = B + 1e-4

        erro_Z = 1.0

        contador_Z = 0


        while erro_Z >= tolerancia_Z and contador_Z < max_iteracoes:

            derivada = dF(Z, A, B)

            Z_novo = Z - F(Z, A, B)/derivada

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

            Z_novo = Z - F(Z, A, B)/derivada

            erro_Z = abs(Z_novo - Z)

            Z = Z_novo

            contador_Z += 1


        Z_vap = Z


        
        # cálculo da diferença de fugacidade
        

        # calcula novamente g no limite inferior do intervalo
        # usando os parâmetros correspondentes a P_inferior

        g_inferior = (
            ln_phi(
                Z_liq,
                a*P_inferior/(R**2*T**2),
                b*P_inferior/(R*T)
            )
            -
            ln_phi(
                Z_vap,
                a*P_inferior/(R**2*T**2),
                b*P_inferior/(R*T)
            )
        )


        # diferença de fugacidade no ponto médio

        g_meio = (
            ln_phi(Z_liq, A, B)
            -
            ln_phi(Z_vap, A, B)
        )


        # verifica em qual metade do intervalo está a solução

        if g_inferior * g_meio < 0:

            P_superior = P_meio

        else:

            P_inferior = P_meio


        # erro associado ao intervalo de pressão

        erro = abs(P_superior - P_inferior)

        contador += 1


    # pressão de saturação calculada pela equação de van der Waals

    P_sat = (P_inferior + P_superior)/2


    # conversão de MPa para kPa

    P_sat_kPa = P_sat*1000


    # valor experimental correspondente à temperatura

    P_nist = nist_data[T_celsius]


    # cálculo do erro percentual em relação ao valor experimental

    erro_percentual = (abs(P_sat_kPa - P_nist)/P_nist)*100


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


# Interlúdio: discussão dos resultados

# Os resultados mostram que a pressão de saturação aumenta com a
# temperatura, acompanhando a tendência observada nos dados do NIST.
# entretanto, os valores calculados pela equação de van der Waals
# apresentam desvios consideráveis em relação aos valores experimentais.
# Acreditamos, sobretudo, que essa diferença releaciona-se às simplificações da equação de
# van der Waals para representar as interações intermoleculares da água.

# Apesar dos desvios quantitativos, o método reproduz o comportamento
# crescente da pressão de saturação com a temperatura.
# Do ponto de vista numérico, o método de Newton-Raphson foi utilizado
# para determinar as raízes da equação cúbica em Z, enquanto a bisseção
# foi utilizada para determinar a pressão de saturação a partir da
# condição de igualdade das fugacidades entre as fases.




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