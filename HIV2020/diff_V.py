import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.integrate import odeint
import networkx as nx
#from NORMtree import N
#from NORMtree import W
#from NORMtree import N_1
N = 7
N_1 = N//3
W = [[0, 0, 1, 1, 0, 0, 0, 0], [0, 0, 1, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 1, 0, 0], [0, 0, 0, 0, 1, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1], [1, 1, 0, 0, 0, 0, 0, 0]]

#1. параметры из статьи
Beta = 2.4 * 10 ** -7
Omega = 0.8 * 10 ** -7
Lambda = 5 * 10 ** 5
d = 0.01
p = 10 ** 3
c = 30
Delta = 0.7
mBL = 1000
pV = 2
pBL = 25
mLB = pBL * mBL
tauj = 1

#из файла с графом достаем N
print()
print('кол-во узлов:', N)
print()
print('матрица связей')
print(W)


#L - число ребер  вместе с кровью
L = 0
for i in range(len(W)):
    for j in range(len(W[i])):
        L += W[i][j]

print('всего ребер', int(L))

#Lj - выходящие ребра из узлов
Lj = []
for row in W:
    Lj.append(int(sum(row)))
print('число выходящих ребер из узла j', Lj)

# mkj = mkk/Lk  - скорость миграции из k в js равна скорости покидания k делить на колво выходящих сосудов

mjj = 1/tauj #скорость покидания узла j

#m - матрица скоростей
m = np.zeros([N+1,N+1])
for i in range(N):
    for j in range(N):
        if W[i][j] != 0:
            m[i][j] = mjj/Lj[i]

for i in range (N-1): #добавили mjj которое вытекает из j
   m[i][i] = mjj

#добавим соединение с кровью для последнего узла и для первых
for i in range(N_1):
    m[N][i] = mBL/Lj[-1]
m[N][N] = mBL
m[N-1][N-1] = mLB
m[N-1][N] = mLB

print('матрица m')
print(m)

mkj = np.zeros([N+1,N+1])
for i in range(N):
    for j in range(N):
        if W[i][j] != 0:
            mkj[i][j] = mjj/Lj[i]

print('mkj=', mkj)

#вектора для клеток начальные концентрации

print('начальные условия для T')
T_0 = []
for i in range(N):
    T_0.append(Lambda/d)
print(T_0)


print('введите начальные условия для I')
I_0 = []
for i in range(N):
    #I_0.append(int(input()))
    I_0.append(0)
print(I_0)

print('введите начальные условия для V')
V_0 = []
for i in range(N):
    #V_0.append(int(input()))
    V_0.append(10)
print(V_0)

#начальные концентрации в крови, пока по нулям
TB_0 = 0
IB_0 = 0
VB_0 = 10


#вектор начальных концентраций, значения для крови в конце
y0 = T_0 + I_0 + V_0
y0.append(TB_0)
y0.append(IB_0)
y0.append(VB_0)
print(y0)

ts = np.linspace(0, 6, 10000)  # time grid

# solve the system dy/dt = f(y, t)
def f(y, t):
    T_I_V = [] # тут T I V плюс TB, IB, VB
    for i in range(N * 3 + 3):
        T_I_V.append(y[i])

    F = []

    for i in range(N):
        F.append(Lambda - (Beta * T_I_V[i + 2 * (N)] + Omega * T_I_V[i + (N)]) * T_I_V[i] - (d + m[i][i]) * T_I_V[i] + (m[N][i]) * T_I_V[3 * (N)])
        for k in range(N):
            F[-1] = F[-1] + mkj[k][i] * T_I_V[k]

    for i in range(N):
        F.append((Beta * T_I_V[i + 2 * (N)] + Omega * T_I_V[i + N]) * T_I_V[i] - (Delta + m[i][i]) * T_I_V[i + N] + (m[N][i]) * T_I_V[3 * (N) + 1 ])
        for k in range(N):
            F[-1] = F[-1] + mkj[k][i] * T_I_V[k + N]

    for i in range(N):
        F.append(p * T_I_V[i + N] - (c + pV * m[i][i]) * T_I_V[i + 2 * (N)] + pV * (m[N][i]) * T_I_V[3 * (N)+2])
        for k in range(N):
            F[-1] = F[-1] + mkj[k][i] * pV * T_I_V[k + 2 * N]

    F.append(- (d + mBL) * T_I_V[3 * (N)] + mLB * T_I_V[N - 1])
    F.append(- (Delta + mBL) * T_I_V[3 * (N) + 1 ] + mLB * T_I_V[(2 * N) - 1])
    F.append(- (c + pV * mBL) * T_I_V[3 * (N) + 2] + pV * mLB * T_I_V[(3 * N) - 1])

    return F


#cамо решение системы

soln = odeint(f, y0, ts)
#for i in range(3*N):
 #   plt.plot(ts, soln[:,i], label='TVI'+str(i), color = (0.3 + i * 0.02 , 0.3, 0.1 + i * 0.02))
#plt.xlabel('Time, days')
#plt.ylabel('T cell [cells/ml], virus [virion/ml]')
#plt.legend(loc=0)

#for i in range(len(t_history_average)):
#    t_history_minimal[i] = np.min(soln[i, 0:N]);
#    v_history_minimal[i] = np.min(soln[i, N:2*N]);
#    i_history_minimal[i] = np.min(soln[i, 2*N:3*N]);
#    t_history_maximal[i] = np.max(soln[i, 0:N]);
#    v_history_maximal[i] = np.max(soln[i, N:2*N]);
#    i_history_maximal[i] = np.max(soln[i, 2*N:3*N]);

#for i in range(N):
#    plt.plot(ts, soln[:,i], label='T'+str(i), color = 'red')
#for i in range(N, 2*N):
#    plt.plot(ts, soln[:, i], label='V' + str(i), color = 'green')
#for i in range(2*N, 3*N):
#    plt.plot(ts, soln[:, i], label='I' + str(i), color = 'blue')


#ТУТ ВЫЧИСЛЯЕМ СРЕДНЕЕ, МИНИМУМ И МАКСИМУМ
t_history_average = np.sum(soln[:, 0:N], axis=1) / N;
i_history_average = np.sum(soln[:, N:2*N], axis=1) / N;
v_history_average = np.sum(soln[:, 2*N:3*N], axis=1) / N;

t_history_minimal = np.min(soln[:, 0:N], axis=1);
i_history_minimal = np.min(soln[:, N:2*N], axis=1);
v_history_minimal = np.min(soln[:, 2*N:3*N], axis=1);

t_history_maximal = np.max(soln[:, 0:N], axis=1);
i_history_maximal = np.max(soln[:, N:2*N], axis=1);
v_history_maximal = np.max(soln[:, 2*N:3*N], axis=1);

plt.plot(ts, np.log10(t_history_minimal), label='T_min', color = 'red', linestyle = 'dashed')
plt.plot(ts, np.log10(t_history_average), label='T_avg', color = 'red')
plt.plot(ts, np.log10(t_history_maximal), label='T_max', color = 'red', linestyle = 'dashed')

plt.plot(ts, np.log10(i_history_minimal), label='I_min', color = 'green', linestyle = 'dashed')
plt.plot(ts, np.log10(i_history_average), label='I_avg', color = 'green')
plt.plot(ts, np.log10(i_history_maximal), label='I_max', color = 'green', linestyle = 'dashed')

plt.plot(ts, np.log10(v_history_minimal), label='V_min', color = 'blue', linestyle = 'dashed')
plt.plot(ts, np.log10(v_history_average), label='V_avg', color = 'blue')
plt.plot(ts, np.log10(v_history_maximal), label='V_max', color = 'blue', linestyle = 'dashed')

plt.title("cреднее")
plt.xlabel('Time, days')
plt.ylabel('T cell [cells/ml], virus [virion/ml] in LS')
plt.legend(loc=0)

plt.show()

plt.figure()
for i in range(N):
    plt.plot(ts, np.log10(soln[:,i]), label='T'+str(i), color = 'red')
for i in range(N, 2*N):
    plt.plot(ts, np.log10(soln[:, i]), label='V' + str(i), color = 'green')
for i in range(2*N, 3*N):
    plt.plot(ts, np.log10(soln[:, i]), label='I' + str(i), color = 'blue')

plt.xlabel('Time, days')
plt.ylabel('T cell [cells/ml], virus [virion/ml] in LS, log')
plt.legend(loc=0)


#в крови
plt.figure()
plt.plot(ts, soln[:,3*N], label='T', color = (0.5 , 0.2, 0.3 ))
plt.plot(ts, soln[:, 3*N+1], label='V' , color=(0.75, 0.75, 0))
plt.plot(ts, soln[:,3*N+2], label='I' , color=(0, 0, 1))

plt.xlabel('Time, days')
plt.ylabel('T cell [cells/ml], virus [virion/ml] in Blood')
plt.legend(loc=0)
plt.show()

#кровь лог
plt.figure()
plt.plot(ts, np.log10(soln[:,3*N]), label='T', color = 'red')
plt.plot(ts, np.log10(soln[:, 3*N+1]), label='I', color = 'green')
plt.plot(ts, np.log10(soln[:,3*N+2]), label='V', color = 'blue')
plt.xlabel('Time')
plt.ylabel('T cell [cells/ml], virus [virion/ml] in Blood, log')
plt.legend(loc=0)
plt.show()


