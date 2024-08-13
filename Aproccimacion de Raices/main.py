def bisectionFor(f, a, b, n):
    if f(a) * f(b) < 0:
        for i in range(n):
            root = (a + b) / 2
            test = f(a) * f(root)
            if f(root) == 0:
                return f'Root found at {root}'
            elif test > 0:
                a = root
            else:
                b = root
        return f'Root approximates to {a}'
    return f'There is no root between {a} and {b}'


def f(x):
    return x ** 3 - 1


print(bisectionFor(f, 1.5, 3, 10))


