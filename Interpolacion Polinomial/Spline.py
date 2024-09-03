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
x = np.array([1, 2, 3])
y = np.array([1, 0.5, 0.3333333])

# # Crear splines
# linear = linear_spline(x, y)
# quadratic = quadratic_spline(x, y)
# cubic = cubic_spline(x, y)
natural = natural_cubic_spline(x, y)

# Crear puntos para graficar
x_plot = np.linspace(0, 5, 200)

# Configurar el gráfico
fig, axs = plt.subplots(2, 2, figsize=(15, 12))
fig.suptitle('Comparación de diferentes tipos de Splines', fontsize=16)

# Función para graficar cada spline
def plot_spline(ax, x_plot, spline, title):
    ax.plot(x, y, 'ro', label='Puntos originales')
    ax.plot(x_plot, spline(x_plot), label='Spline')
    ax.set_title(title)
    ax.legend()
    ax.grid(True)

# Graficar cada spline
# plot_spline(axs[0, 0], x_plot, linear, 'Spline Lineal')
# plot_spline(axs[0, 1], x_plot, quadratic, 'Spline Cuadrático')
# plot_spline(axs[1, 0], x_plot, cubic, 'Spline Cúbico')
plot_spline(axs[1, 1], x_plot, natural, 'Spline Cúbico Natural')

plt.tight_layout()
plt.show()

# Evaluar cada spline en un punto de ejemplo
x_eval = 2.5
print(f"Valores interpolados en x = {x_eval}:")
# print(f"Spline Lineal: {linear(x_eval):.4f}")
# print(f"Spline Cuadrático: {quadratic(x_eval):.4f}")
# print(f"Spline Cúbico: {cubic(x_eval):.4f}")
print(f"Spline Cúbico Natural: {natural(x_eval):.4f}")