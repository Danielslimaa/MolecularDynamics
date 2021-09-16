# AUTOR: DANIEL SOUZA LIMA 15 de setembro de 2021
#tradução do código em .c feito pelo Prof. Dr. Lucas Nicolao

#/* **************************************************************************************
#Algortimo básico simulação dinâmica molecular no ensemble NVE.
#Potencial de pares de Lennard-Jones em D=3, levando em conta todos pares.

# Parâmetros da simulação: temperatura, densidade (rho), passo de tempo
#(dt), número de partículas (N), dimensão do sistema (D)

# Descrições das funções:
#1) force(double r[][D], double a[][D])
#- calcula a força resultante em cada partícula/direção cartesiana, e armazena em 'a'
#2) vverlet(double r[][D], double v[][D], double a[][D])
#- atualização das posições 'r' e velocidades 'v' de acordo com velocity Verlet
#- devolve o valor da energia potencial por partícula
#3) measures(double r[][D], double v[][D], double *energia, double *temp, double *pressao)
#- mede energia total por partícula, temperatura cinética e pressão virial
#- lembrando que 3*temp/2 = energia cinética e energia - 3*temp/2 = energia potencial
#4) double energiacin(double v[][D])
#- devolve energia cinética por partícula
#5) overrelax(double r[][D], double a[][D])
#- suaviza forças de uma condição inicial aleatórias segundo dr/dt = -grad U
#6) initial3D(double r[][D], double v[][D], int qual)
#- condições iniciais aleatórias (qual=0) e cristal cúbico (qual=1) para posições
#- velocidades aleatórias de acordo com parâmetro temperatura
#7) reescalavT(double v[][D], double Talvo)
#- rescala velocidades para atingir temperatura alvo Talvo
#8) reescalarRho(double r[][D], double rho_alvo)
#- reescala posições, tamanho da caixa e densidade para mudar densidade para rho_alvo
#9) printXYZ(double r[][D])
#- imprime (na tela) configurações p/ compor arquivo xyz. 1a partícula cor diferente.
#************************************************************************************** */


import numpy as np
import math
import random

global temperatura, rho, dt, N, D, L, V

temperatura = (1.0)
rho = (0.80)

dt = (0.001)

N = int(100)
D = int(3)

L = pow((N)/rho, 1.0/(D))
V = pow(L, (D))

print("Hello, world!")


def force(r, a):

    a = np.zeros((N, D))
    dr = np.zeros(D)
    en = 0

    for i in range(0, N-1):
        for j in range(i + 1, N):

            d2 = 0
            for n in range(0, D):
                dr[n] = r[i, n] - r[j, n]
                dr[n] = dr[n] - L * math.floor(dr[n] / L + 0.5)
                d2 += pow(dr[n], 2)

            r2 = 1.0 / d2
            r6 = pow(r2, 3)
            ff = 48.0 * r2 * r6 * (r6 - 0.5)

            for n in range(0, D):
                a[i, n] += ff*dr[n]
                a[j, n] -= ff*dr[n]
            en += 4.0 * r6 * (r6 - 1.0)

    return a, en / N


def vverlet(r, v, a):

    energia = 0
    a, energia = force(r, a)
    for i in range(0, N):
        for n in range(0, D):
            v[i, n] += 0.5*dt*a[i, n]
            r[i, n] += dt*v[i, n]

    a, energia = force(r, a)
    for i in range(0, N):
        for n in range(0, D):
            v[i, n] += 0.5*dt*a[i, n]

    return r, v, a


def measures(r, v, energia, temp, pressao):
    dr = np.zeros(D)
    d2 = 0
    r2 = 0
    r6 = 0
    ff = 0
    virial = 0.
    sumv2 = 0.
    en = 0

    for i in range(0, N - 1):

        for n in range(0, D):
            sumv2 += v[i, n] * v[i, n]

        for j in range(i + 1, N):
            d2 = 0
            for n in range(0, D):
                dr[n] = r[i, n] - r[j, n]
                dr[n] = dr[n] - L * math.floor((dr[n] / L) + 0.5)
                d2 += pow(dr[n], 2)
            r2 = 1.0 / d2
            r6 = pow(r2, 3)
            ff = 48.0 * r2 * r6 * (r6 - 0.5)

            for n in range(0, D):
                virial += ff * dr[n] * dr[n]
            en += 4.0 * r6 * (r6 - 1.0)

    for n in range(0, D):
        sumv2 += v[N - 1, n] * v[N-1, n]

    energia = sumv2 / (2.0 * (N)) + (en / (N))
    temp = sumv2 / ((D) * (N))
    pressao = temp * rho + virial / ((V) * (D))

    return r, v, energia, temp, pressao


def energiacin(v):
    K = 0.0

    for i in range(0, N):
        for n in range(0, D):
            K += v[i, n] * v[i, n]

    return K / (2.0 * N)


def overrelax(r, a):

    Dt = 0.1
    energia = 0.0

    a, energia = force(r, a)
    for i in range(0, N):
        norma = 0.0
        for n in range(0, D):
            norma += pow(a[i, n], 2)
        norma = math.sqrt(norma)
        for n in range(0, D):
            r[i, n] += Dt * a[i, n] / norma

    return


def initial3D(r, v, qual):

    somav2 = 0
    somav = np.zeros(D)

    if qual == 1:
        Nsites = int(round(pow(N, 1.0 / D)))
        dx = L / (Nsites)
        for i in range(0, Nsites):
            for j in range(0, Nsites):
                for k in range(0, Nsites):
                    ii = k + Nsites * (j + i * Nsites)
                    if ii < N:
                        r[ii, 0] = (i + 0.5) * dx
                        r[ii, 1] = (j + 0.5) * dx
                        r[ii, 2] = (k + 0.5) * dx
    else:
        for i in range(0, N):
            for n in range(0, D):
                r[i, n] = random.random() * L

    for i in range(0, N):
        for n in range(0, D):
            v[i, n] = random.random() - 0.5
            somav[n] += v[i, n]
            somav2 += math.pow(v[i, n], 2)

    for n in range(0, D):
        somav[n] /= N
    somav2 /= N
    fac = math.sqrt(D * temperatura / somav2)
    for i in range(0, N):
        for n in range(0, D):
            v[i, n] = (v[i, n] - somav[n]) * fac

    return v


def reescalavT(v, Talvo):
    temp, somav2 = 0

    for i in range(0, N):
        for n in range(0, D):
            somav2 += pow(v[i, n], 2)

    somav2 /= N
    temp = somav2 / D
    fac = math.sqrt(Talvo / temp)
    for i in range(0, N):
        for n in range(0, D):
            v[i, n] *= fac

    return v


def reescalarRho(r, rho_alvo):

    fac = pow(rho / rho_alvo, 1.0 / D)

    for i in range(0, N):
        for n in range(0, D):
            r[i, n] *= fac
    rho = rho_alvo

    return rho


def printXYZ(r):

    print(str(N) + "\n\n")
    for i in range(0, N):
        if i == 0:
            aux = 'B'
        else:
            aux = 'A'
        print(str(aux))
        for n in range(0, D):
            print(str(r[i][n] - L * math.floor(r[i][n]/L + 0.5) + L/2.) + "\n")


#Impressora para o OVITO
def impressora_video(r, N, tempo):
    
    buffer = "ITEM: TIMESTEP\n" + str(tempo) + "\n" + "ITEM: NUMBER OF ATOMS\n" + str(N) + "\n" + "ITEM: BOX BOUNDS ss ss pp\n" 
    buffer += "-1" + " " + "1\n" + "-1" + " " +  "1\n"
    buffer += "-1" + " " + "1\n" + "ITEM: ATOMS id x y z" + "\n"

    with open("particles.dump", 'a') as file_object:
        file_object.write(buffer)
        for i in range(N):
            buffer2 = str(int(i)) + "\t" +  str(round(r[i,0],3)) + "\t" + str(round(r[i,1],3)) + "\t" + str(round(r[i,2],3)) + "\n"
            file_object.write(buffer2)


#Lê as condições iniciais --------------- COMENTE A CHAMADA da FUNÇÃO "initi3D" SE FOR USÁ-lA
def ler_CondicoesIniciais(nome_do_arquivo):

    with open(nome_do_arquivo, "r") as file_object:
        dados = np.loadtxt(file_object)

    N = dados.shape[0]
    r = np.zeros((dados.shape[0],int(dados.shape[1]/2)))
    v = np.zeros((dados.shape[0],int(dados.shape[1]/2)))

    r = dados[:,0:3]
    v = dados[:,3:6]
    
    return r, v, N

# Começo do programa

r = np.zeros((N, D))
v = np.zeros((N, D))
a = np.zeros((N, D))

K = 0.0
U = 0.0
E = 0.0
T = 0.0
P = 0.0

#Condição inicial aleatória
initial3D(r, v, 0) ### COMENTE ISSO SE FOR USAR "ler_CondicoesIniciais("condicoes_iniciais.dat")"
K = energiacin(v)

#Condição inicial espicificada a partir do arquivo "condicoes_iniciais.dat"
#no formato x y z vx vy vz em que cada linha coresponde a uma partícula
#r, v, N = ler_CondicoesIniciais("condicoes_iniciais.dat")

for t in range(0, 5 * N):
    overrelax(r, a)
    a, U = force(r, a)
    if U <= 0:
        break

for t in range(0, 1000):
    r, v, a = vverlet(r, v, a)

    r, v, E, T, P = measures(r, v, E, T, P)

    K = 3.0 * T / 2.0
    U = E - K
    print(str(t) + " " + str(round(K, 6)) + " " + str(round(U, 6)) + " " + str(round(E, 6)) + " " + str(round(T, 6)) + " " + str(round(P, 6)))
    impressora_video(r, N, t)



