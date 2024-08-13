# 1. y=10*cos(x) entre 0 y 2 pi
import numpy as np
from matplotlib import pyplot as plt

x = np.linspace(0, 2*np.pi, 100)
y = 10*np.cos(x)

plt.axvline(x=0, color='r')
plt.axhline(y=0, color='r')
plt.plot(x,y, color='blue')
plt.show()

# 2. z=t^2+1/3 * t entre -10 y 20 con etiquetas de aceleracion y tiempo

t = np.linspace(-10, 20, 100)
z = t**2 + (1/3) * t

plt.axvline(x=0, color='r')
plt.axhline(y=0, color='r')
plt.xlabel('Tiempo')
plt.ylabel('Aceleracion')
plt.plot(t,z, color='blue')
plt.show()

# 3. curva parametrica

t = np.linspace(0, 2*np.pi, 100)
x = np.cos(t)
y = np.sin(t)

plt.plot(x, y)
plt.title('Curva Paramétrica: Círculo')
plt.axis('equal')
plt.axvline(x=0, color='r')
plt.axhline(y=0, color='r')
plt.show()