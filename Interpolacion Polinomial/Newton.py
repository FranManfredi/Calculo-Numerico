import numpy as np

def newton_interpolation_polynomial(x, y):
    n = len(x)
    
    # Calcula la tabla de diferencias divididas
    coef = np.zeros((n, n))
    coef[:,0] = y
    
    for j in range(1, n):
        for i in range(n - j):
            coef[i][j] = (coef[i+1][j-1] - coef[i][j-1]) / (x[i+j] - x[i])
    
    # Obtiene los coeficientes del polinomio
    a = coef[0]
    
    # Crea el polinomio
    def p(x_val):
        n = len(a) - 1
        temp = a[n]
        for k in range(n - 1, -1, -1):
            temp = temp * (x_val - x[k]) + a[k]
        return temp
    
    # Crea una representación en cadena del polinomio
    def poly_string():
        terms = []
        for i in range(n):
            if a[i] != 0:
                term = f"{a[i]:.4f}"
                for j in range(i):
                    term += f"*(x - {x[j]})"
                terms.append(term)
        return " + ".join(terms)
    
    p.coef = a
    p.string = poly_string
    
    return p

# Ejemplo de uso
x = [0, 1, 2, 3]
y = [1, 2, 4, 8]

p = newton_interpolation_polynomial(x, y)

print("Polinomio interpolador:")
print(p.string())

# Evaluar el polinomio en un punto
x_interp = 1.5
resultado = p(x_interp)
print(f"\nEl valor interpolado en x = {x_interp} es {resultado}")
