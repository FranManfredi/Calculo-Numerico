import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, CubicSpline

def linear_spline(x, y):
    return interp1d(x, y, kind='linear')

def quadratic_spline(x, y):
    return interp1d(x, y, kind='quadratic')

def cubic_spline(x, y):
    return interp1d(x, y, kind='cubic')

def natural_cubic_spline(x, y):
    return CubicSpline(x, y, bc_type='natural')

# Datos de ejemplo
x = np.array([0, 1, 2, 3, 4, 5])
y = np.array([2, 1, 4, 3, 5, 4])

# Crear splines
linear = linear_spline(x, y)
quadratic = quadratic_spline(x, y)
cubic = cubic_spline(x, y)
natural = natural_cubic_spline(x, y)

# Crear puntos para graficar
x_plot = np.linspace(0, 5, 200)

# Configurar el gráfico
plt.figure(figsize=(12, 8))
plt.title('Comparación de diferentes tipos de Splines', fontsize=16)

# Graficar cada spline
plt.plot(x, y, 'ko', markersize=10, label='Puntos originales')
plt.plot(x_plot, linear(x_plot), 'r-', label='Spline Lineal')
plt.plot(x_plot, quadratic(x_plot), 'g-', label='Spline Cuadrático')
plt.plot(x_plot, cubic(x_plot), 'b-', label='Spline Cúbico')
plt.plot(x_plot, natural(x_plot), 'm-', label='Spline Cúbico Natural')

plt.xlabel('x', fontsize=12)
plt.ylabel('y', fontsize=12)
plt.legend(fontsize=10)
plt.grid(True)
plt.tight_layout()
plt.show()

# Evaluar cada spline en un punto de ejemplo
x_eval = 2.5
print(f"Valores interpolados en x = {x_eval}:")
print(f"Spline Lineal: {linear(x_eval):.4f}")
print(f"Spline Cuadrático: {quadratic(x_eval):.4f}")
print(f"Spline Cúbico: {cubic(x_eval):.4f}")
print(f"Spline Cúbico Natural: {natural(x_eval):.4f}")