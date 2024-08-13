import numpy as np
import matplotlib.pyplot as plt

x = np.arange(-0.5, 1.5, 0.0000001)
y = x * np.cos(x)
# x1 =
y1 = x * np.sin(x)

plt.plot(x, y, color='blue')
plt.plot(x, y1, color='red')
plt.axvline(x=0, color='black')
plt.axhline(y=0, color='black')
plt.show()

x = np.linspace(0, 30, 100)  # el num te da el valor de secciones que te da la funcion
y = np.exp(-0.1 * x) * np.cos(x)
plt.plot(x, y, "r--", color='blue')
plt.show()

