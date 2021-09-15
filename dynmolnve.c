/* **************************************************************************************
Algortimo básico simulação dinâmica molecular no ensemble NVE.
Potencial de pares de Lennard-Jones em D=3, levando em conta todos pares.

# Parâmetros da simulação: temperatura, densidade (rho), passo de tempo
(dt), número de partículas (N), dimensão do sistema (D)

# Descrições das funções:
1) force(double r[][D], double a[][D])
- calcula a força resultante em cada partícula/direção cartesiana, e armazena em 'a'
2) vverlet(double r[][D], double v[][D], double a[][D])
- atualização das posições 'r' e velocidades 'v' de acordo com velocity Verlet
- devolve o valor da energia potencial por partícula
3) measures(double r[][D], double v[][D], double *energia, double *temp, double *pressao)
- mede energia total por partícula, temperatura cinética e pressão virial
- lembrando que 3*temp/2 = energia cinética e energia - 3*temp/2 = energia potencial
4) double energiacin(double v[][D])
- devolve energia cinética por partícula
5) overrelax(double r[][D], double a[][D])
- suaviza forças de uma condição inicial aleatórias segundo dr/dt = -grad U
6) initial3D(double r[][D], double v[][D], int qual)
- condições iniciais aleatórias (qual=0) e cristal cúbico (qual=1) para posições
- velocidades aleatórias de acordo com parâmetro temperatura
7) reescalavT(double v[][D], double Talvo)
- rescala velocidades para atingir temperatura alvo Talvo
8) reescalarRho(double r[][D], double rho_alvo)
- reescala posições, tamanho da caixa e densidade para mudar densidade para rho_alvo
9) printXYZ(double r[][D])
- imprime (na tela) configurações p/ compor arquivo xyz. 1a partícula cor diferente.
************************************************************************************** */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define temperatura 1.0
double rho = 0.80;

#define dt 0.001

#define N 100
#define D 3

#define L (pow(N / rho, 1.0 / D))
#define V (pow(L, D))

#define FRANDOM (random() / (RAND_MAX + 1.0))

double force(double r[][D], double a[][D])
{
  int i, j, n;
  double dr[D], d2, r2, r6, ff, en = 0;

  for (i = 0; i < N; i++)
    for (n = 0; n < D; n++)
      a[i][n] = 0;

  for (i = 0; i < N - 1; i++)
  {
    for (j = i + 1; j < N; j++)
    {

      d2 = 0;
      for (n = 0; n < D; n++)
      {
        dr[n] = r[i][n] - r[j][n];
        dr[n] = dr[n] - L * floor(dr[n] / L + 0.5);
        d2 += pow(dr[n], 2);
      }
      r2 = 1 / d2;
      r6 = r2 * r2 * r2;
      ff = 48 * r2 * r6 * (r6 - 0.5);

      for (n = 0; n < D; n++)
      {
        a[i][n] += ff * dr[n];
        a[j][n] -= ff * dr[n];
      }
      en += 4 * r6 * (r6 - 1);
    }
  }

  return en / N; 
}

void vverlet(double r[][D], double v[][D], double a[][D])
{
  int i, n;

  force(r, a);
  for (i = 0; i < N; i++)
  {
    for (n = 0; n < D; n++)
    {
      v[i][n] += 0.5 * dt * a[i][n];
      r[i][n] += dt * v[i][n];
    }
  }
  force(r, a);
  for (i = 0; i < N; i++)
    for (n = 0; n < D; n++)
      v[i][n] += 0.5 * dt * a[i][n];

  return;
}

void measures(double r[][D], double v[][D], double *energia, double *temp, double *pressao)
{
  int i, j, n;
  double dr[D], d2, r2, r6, ff;
  double virial = 0, sumv2 = 0, en = 0;

  for (i = 0; i < N - 1; i++)
  {

    for (n = 0; n < D; n++)
      sumv2 += v[i][n] * v[i][n];

    for (j = i + 1; j < N; j++)
    {
      d2 = 0;
      for (n = 0; n < D; n++)
      {
        dr[n] = r[i][n] - r[j][n];
        dr[n] = dr[n] - L * floor(dr[n] / L + 0.5);
        d2 += pow(dr[n], 2);
      }
      r2 = 1 / d2;
      r6 = r2 * r2 * r2;
      ff = 48 * r2 * r6 * (r6 - 0.5);

      for (n = 0; n < D; n++)
        virial += ff * dr[n] * dr[n];
      en += 4 * r6 * (r6 - 1);
    }
  }

  for (n = 0; n < D; n++)
    sumv2 += v[N - 1][n] * v[N - 1][n];

  *energia = sumv2 / 2 / N + en / N;
  *temp = sumv2 / D / N;
  *pressao = *temp * rho + virial / V / D;

  return;
}

double energiacin(double v[][D])
{
  int i, n;
  double K = 0;

  for (i = 0; i < N; i++)
    for (n = 0; n < D; n++)
      K += v[i][n] * v[i][n];

  return K / 2.0 / N;
}

void overrelax(double r[][D], double a[][D])
{
  int i, n;
  double norma, Dt = 0.1;
  /* dr/dt = -grad U */

  force(r, a);
  for (i = 0; i < N; i++)
  {
    norma = 0;
    for (n = 0; n < D; n++)
      norma += pow(a[i][n], 2);
    norma = sqrt(norma);
    for (n = 0; n < D; n++)
      r[i][n] += Dt * a[i][n] / norma;
  }

  return;
}

void initial3D(double r[][D], double v[][D], int qual)
{
  int i, j, k, n, ii;
  int Nsites;
  double dx, somav2 = 0, somav[D] = {0}, fac;

  if (qual == 1)
  {
    Nsites = (int)round(pow(N, 1.0 / D));
    dx = L / (double)Nsites;
    for (i = 0; i < Nsites; i++)
    {
      for (j = 0; j < Nsites; j++)
      {
        for (k = 0; k < Nsites; k++)
        {
          ii = k + Nsites * (j + i * Nsites);
          if (ii < N)
          {
            r[ii][0] = (i + .5) * dx;
            r[ii][1] = (j + .5) * dx;
            r[ii][2] = (k + .5) * dx;
          }
        }
      }
    }
  }
  else
  {
    for (i = 0; i < N; i++)
      for (n = 0; n < D; n++)
        r[i][n] = FRANDOM * L;
  }

  for (i = 0; i < N; i++)
  {
    for (n = 0; n < D; n++)
    {
      v[i][n] = FRANDOM - 0.5;
      somav[n] += v[i][n];
      somav2 += pow(v[i][n], 2);
    }
  }
  for (n = 0; n < D; n++)
    somav[n] /= N;
  somav2 /= N;
  fac = sqrt(D * temperatura / somav2);
  for (i = 0; i < N; i++)
    for (n = 0; n < D; n++)
      v[i][n] = (v[i][n] - somav[n]) * fac;

  return;
}

void reescalavT(double v[][D], double Talvo)
{
  int i, n;
  double temp, somav2 = 0, fac;

  for (i = 0; i < N; i++)
    for (n = 0; n < D; n++)
      somav2 += pow(v[i][n], 2);

  somav2 /= N;
  temp = somav2 / D;
  fac = sqrt(Talvo / temp);
  for (i = 0; i < N; i++)
    for (n = 0; n < D; n++)
      v[i][n] *= fac;

  return;
}

void reescalarRho(double r[][D], double rho_alvo)
{
  int i, n;
  double fac;

  fac = pow(rho / rho_alvo, 1.0 / D);

  for (i = 0; i < N; i++)
    for (n = 0; n < D; n++)
      r[i][n] *= fac;

  rho = rho_alvo;

  return;
}

void printXYZ(double r[][D])
{
  int i, n;

  printf("%d\n\n", N);
  for (i = 0; i < N; i++)
  {
    printf("%c ", (i == 0 ? 'B' : 'A'));
    for (n = 0; n < D; n++)
      printf("%f ", r[i][n] - L * floor(r[i][n] / L + 0.5) + L / 2.);
    printf("\n");
  }

  return;
}

int main()
{
  int t;
  double r[N][D], v[N][D], a[N][D];
  double K, U, E, T, P;

  //exemplo condição inicial aleatória
  initial3D(r, v, 0);
  K = energiacin(v);
  for (t = 0; t < 5 * N; t++)
  {
    overrelax(r, a);
    U = force(r, a);
    if (U <= 0)
      break;
    //printf("%d %f %f %f\n",t,U,K,U+K);
  }

  for (t = 0; t < 1000; t++)
  {
    vverlet(r, v, a);

    measures(r, v, &E, &T, &P);
    K = 3 * T / 2;
    U = E - K;
    printf("%d %f %f %f %f %f\n", t, K, U, E, T, P);
  }

  return 0;
}