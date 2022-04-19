import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-np.pi, np.pi, 100)
#y = np.power(x, 2) - np.log(x) - 2 * np.cos(x) - 1
y = np.tan(0.4*x+0.4)-np.power(x, 2)

fig, ax = plt.subplots()

ax.plot(x, y)
ax.grid()

fig.set_figwidth(8)
fig.set_figheight(8)
plt.show()

# dyhotomy
def fun(x):
    #return np.tan(0.4*x+0.4)-np.power(x, 2)
    return (1/x) - 7
    
def half_del(a, b, eps):
    i = 0
    if fun(a)*fun(b) < 0:
        while abs(a-b) > eps:
            i += 1
            c = (a+b)/2
            if fun(a)*fun(c) < 0:
                b = c
            else:
                a = c
    else:
        print('Access denied')
        return None
    return a, i, eps

a, b = float(input('Начало отрезка: ')), float(input('Конец отрезка: '))
root, iterations, eps = half_del(a, b, 10 ** int(input('Степень точности = ')))

print()
print('Уравнение: tg(0.4x+0.4)-x^2=0')
print('Отрезок:', '[', a, ',', b, ']')
print('Точность:', eps)
print('Метод половинного деления')
print('x =', root)
print('n =', iterations)
print('Проверка:', fun(root))

# newton
import math
def fun(x):
    return math.tan(0.4*x+0.4)-x**2

def dfdx(x):
    return -2*x+0.4*(1/(math.cos(0.4*x+0.4))**2)

def newton_method(a, b, eps):
    i = 0
    c = a-(fun(a)*(b-a))/(fun(b)-fun(a))
    if fun(a)*fun(b) < 0: # формируем начальное приближение, так как не задана вторая производная функции
        i_x = [a]
    elif fun(a)*fun(b) > 0:
        i_x = [b]
    elif fun(c) == 0:
        i_x = [c]
        
    while a > eps:
        i += 1
        i_x.append(i_x[i-1]-(fun(i_x[i-1])/dfdx(i_x[i-1])))
        a = abs(i_x[i-1]-i_x[i])
    return i_x[-1], i, eps

a, b = float(input('Начало отрезка: ')), float(input('Конец отрезка: '))
root, iterations, eps = newton_method(a, b, 10 ** int(input('Степень точности = ')))

print()
print('Уравнение: tg(0.4x+0.4)-x^2=0')
print('Отрезок:', '[', a, ',', b, ']')
print('Точность:', eps)
print('Метод касательных')
print('x =', root)
print('n =', iterations)
print('Проверка:', fun(root))

# secant method
import numpy as np
def fun(x):
    return np.tan(0.4*x+0.4)-np.power(x, 2)
    
def secant(a, b, eps):
    i = 0
    e = abs(a-b)
    i_x = [0]
    if fun(a)*fun(b) < 0:
        while e > eps:
            i += 1
            c = a-(fun(a)*(b-a))/(fun(b)-fun(a))
            if fun(a)*fun(c) < 0:
                b = c
            else:
                a = c
            i_x.append(c)
            e = abs(i_x[i]-i_x[i-1])
    else:
        print('Access denied')
        return None
    return a, i, eps

a, b = float(input('Начало отрезка: ')), float(input('Конец отрезка: '))
root, iterations, eps = secant(a, b, 10 ** int(input('Степень точности = ')))

print()
print('Уравнение: tg(0.4x+0.4)-x^2=0')
print('Отрезок:', '[', a, ',', b, ']')
print('Точность:', eps)
print('Метод хорд')
print('x =', root)
print('n =', iterations)
print('Проверка:', fun(root))

# secant method
import math
def fun(x):
    return math.tan(0.4*x+0.4)-x**2
    
def second_diff(x):
    return 0.32*(1/(math.cos(0.4*x+0.4))**2)*math.tan(0.4*x+0.4)-2
    
def method_chord(a, b, eps):
    if fun(a)*fun(b) < 0:
        if fun(a)*second_diff(a) > 0:
            c = a
        if fun(b)*second_diff(b) > 0:
            c = b
        if fun(a)*second_diff(a) < 0:
            x = a
        if fun(b)*second_diff(b) < 0:
            x = b
        i = 1
        delta = fun(x) * (x - c) / (fun(x) - fun(c))
        x = x - delta
        while abs(delta) > eps:
            delta = fun(x) * (x - c) / (fun(x) - fun(c))
            x = x - delta
            i += 1
        return x, i, eps
    else:
        print('Access denied')
        return None

a, b = float(input('Начало отрезка: ')), float(input('Конец отрезка: '))
root, iterations, eps = method_chord(a, b, 10 ** int(input('Степень точности = ')))

print()
print('Уравнение: tg(0.4x+0.4)-x^2=0')
print('Отрезок:', '[', a, ',', b, ']')
print('Точность:', eps)
print('Метод хорд')
print('x =', root)
print('n =', iterations)
print('Проверка:', fun(root))

# fixed point iteration
import math
def fun(x):
    return math.tan(0.4*x+0.4)-x**2

def equal_fun(x):
    return math.sqrt(math.tan(0.4*x+0.4))

def fixed_point(a, b, eps):
    if fun(a)*fun(b) < 0:
        i = 0
        x0 = float(input('Начальное приближение: '))
        if a <= x0 <= b:
            while True:
                i += 1
                x1 = equal_fun(x0)
                if math.fabs(x1 - x0) < eps: return x1, i, eps
                x0 = x1
        else:
            print('Начальное приближение неудовлетворяет отрезку')
            return None
    else:
        print('Access denied')
        return None

a, b = float(input('Начало отрезка: ')), float(input('Конец отрезка: '))
root, iterations, eps = fixed_point(a, b, 10 ** int(input('Степень точности = ')))

print()
print('Уравнение: tg(0.4x+0.4)-x^2=0')
print('Отрезок:', '[', a, ',', b, ']')
print('Точность:', eps)
print('Метод последовательных приближений')
print('x =', root)
print('n =', iterations)
print('Проверка:', fun(root))

# Chebyshev method
import math
def fun(x):
    return math.tan(0.4*x+0.4)-x**2

def first_diff(x):
    return -2*x+0.4*(1/(math.cos(0.4*x+0.4))**2)
    
def second_diff(x):
    return 0.32*(1/(math.cos(0.4*x+0.4))**2)*math.tan(0.4*x+0.4)-2

def chebyshev_method(a, b, eps):
    if fun(a)*fun(b) < 0:
        i = 0
        x0 = float(input('Начальное приближение: '))
        if a <= x0 <= b:       
            while True:
                i += 1
                f_n = fun(x0)
                f_1 = first_diff(x0)
                f_2 = second_diff(x0)
                x_n = x0 - (f_n/f_1) - (f_2*f_n**2)/((f_1**3)*2)
                if math.fabs(x_n - x0) < eps: return x_n, i, eps
                x0 = x_n
        else:
            print('Начальное приближение неудовлетворяет отрезку')
            return None
    else:
        print('Access denied')
        return None

a, b = float(input('Начало отрезка: ')), float(input('Конец отрезка: '))
root, iterations, eps = chebyshev_method(a, b, 10 ** int(input('Степень точности = ')))

print()
print('Уравнение: tg(0.4x+0.4)-x^2=0')
print('Отрезок:', '[', a, ',', b, ']')
print('Точность:', eps)
print('Метод Чебышева')
print('x =', root)
print('n =', iterations)
print('Проверка:', fun(root))
