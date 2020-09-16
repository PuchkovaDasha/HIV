import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.integrate import odeint

Beta = 2 * 10 ** (-7)
Omega = 0.8 * 10 ** -7
Lambda = 10 ** 5
d = 0.1
p = 100
c = 5
Delta = 0.5

# initial concentration (начальная концентрация)
T_0 = Lambda/d
I_0 = 0
V_0 = 10

y0 = [T_0, I_0, V_0]

ts2 = np.linspace(0, 30, 1000000)

# solve the system dy/dt = f(y, t)
def f(y, t):
    T = y[0]
    I = y[1]
    V = y[2]

    # the model equations
    f1 = Lambda - Beta * T * V  - d * T
    f2 = Beta * T * V  - Delta * I
    f3 = p * I - c * V

    return [f1, f2, f3]

# тут 0 инфицированных и 10 вирионов
soln = odeint(f, y0, ts2)
TR = soln[:, 0]
IR = soln[:, 1]
VR = soln[:, 2]

plt.figure()
plt.semilogy(ts2, TR, label='T', color = 'red')
plt.semilogy(ts2, IR, label='I', color = 'green')
plt.semilogy(ts2, VR, label='V', color = 'blue')
plt.xlabel('Время (дни)')
plt.ylabel('T (клетки/мл), I (клетки/мл), V (вирионы/мл)')
plt.legend(loc=0)
plt.show()

plt.figure()
plt.subplot(3, 1, 1)
plt.semilogy(ts2, TR, label='T', color = 'red')            # построение графика
plt.ylabel("y1", fontsize=14) # ось ординат
#plt.xlabel('Время (день)')
plt.ylabel('T (клетки/мл)', fontsize=12)

#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.subplot(3, 1, 2)
plt.semilogy(ts2, IR, label='I', color = 'green')               # построение графика
#plt.xlabel('Время (день)')
plt.ylabel('I (клетки/мл)', fontsize=12)
plt.ylim(10**(2), 10**(6))
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.subplot(3, 1, 3)
plt.semilogy(ts2, VR, label='V', color = 'blue')           # построение графика
plt.xlabel('Время (дни)', fontsize=12)
plt.ylabel('V (вирионы/мл)', fontsize=12)
plt.ylim(10**(3), 10**(8))
plt.show()

# тут 0 инфицированных и 10 вирионов
plt.figure()
#plt.plot(ts2, np.log10(TR), label='T', color = 'red')
#plt.plot(ts2, np.log10(IR), label='I', color = 'green')
#plt.plot(ts2, np.log10(VR), label='V', color = 'blue')
#plt.xlabel('Time, days')
#plt.ylabel('T , I, V, cells/ml (log)')
#plt.legend(loc=0)