## AUTOR: DANIEL SOUZA LIMA de setembro de 2021 

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



#### Começo do programa

r = np.zeros((N,D))
v = np.zeros((N,D))
a = np.zeros((N,D))

K = 0.0; U = 0.0; E = 0.0; T = 0.0; P = 0.0

initial3D(r, v, 0)
K = energiacin(v)

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
    print(str(t) + " " + str(round(K,6)) + " " + str(round(U,6)) + " " +
     str(round(E,6)) +  " " + str(round(T,6)) +  " " + str(round(P,6)))

