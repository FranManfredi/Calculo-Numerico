import numpy as np
import matplotlib.pyplot as plt

def lagrange_interpolation(x, y):
    n = len(x)
    
    def L(k, x_val):
        """Función base de Lagrange"""
        return np.prod([(x_val - x[j])/(x[k] - x[j]) for j in range(n) if j != k])
    
    def P(x_val):
        """Polinomio interpolador de Lagrange"""
        return sum(y[k] * L(k, x_val) for k in range(n))
    
    def poly_string():
        """Representación en cadena del polinomio"""
        terms = []
        for k in range(n):
            term = f"{y[k]}"
            for j in range(n):
                if j != k:
                    term += f" * (x - {x[j]}) / ({x[k] - x[j]})"
            terms.append(f"({term})")
        return " + ".join(terms)
    
    P.string = poly_string
    return P

# Ejemplo de uso
x = [-1, 0, 1, 2]
y = [1, -1, 1, 5]

# Crear el polinomio interpolador
P = lagrange_interpolation(x, y)

# Mostrar el polinomio
print("Polinomio interpolador de Lagrange:")
print(P.string())

# Evaluar el polinomio en un punto
x_interp = 0.5
y_interp = P(x_interp)
print(f"\nValor interpolado en x = {x_interp}: {y_interp}")

# Visualización
x_plot = np.linspace(min(x) - 0.5, max(x) + 0.5, 100)
y_plot = [P(xi) for xi in x_plot]

plt.figure(figsize=(10, 6))
plt.plot(x_plot, y_plot, 'b-', label='Polinomio interpolador')
plt.plot(x, y, 'ro', label='Puntos dados')
plt.plot(x_interp, y_interp, 'g*', markersize=10, label='Punto interpolado')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Interpolación de Lagrange')
plt.legend()
plt.grid(True)
plt.show()