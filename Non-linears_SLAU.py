import math
import numpy as np
eps = 10**-6
x0_1=1
x0_2=2

while True:
    system_eq = np.array([[(x0_1**2)-x0_2-1],[x0_1-(x0_2**2)+1]], float)
    #print(system_eq)
    
    jacobian = np.array([[2*x0_1, -1],[1, -2*x0_2]], float)
    jacobian_inv = np.linalg.inv(jacobian)
    
    temp = jacobian_inv.dot(system_eq)
    
    temp = np.array([[x0_1], [x0_2]], float) - temp
    x0_1 = temp[0][0]
    x0_2 = temp[1][0]
    if abs(system_eq[0][0]) <= eps >= abs(system_eq[1][0]):
        print(x0_1, x0_2)
        break
        

eps = 10**-4
x0_1=3.8
x0_2=2
i = 0
while True:
    i += 1
    system_eq = np.array([[x0_1**3+x0_2**3-8*x0_2*x0_1],[x0_1*np.log(x0_2)-x0_2*np.log(x0_1)]], float)
    
    jacobian = np.array([[3*(x0_1**2)-8*x0_2, 3*(x0_2**2)-8*x0_1],[np.log(x0_2)-(x0_2/x0_1), (x0_1/x0_2)-np.log(x0_1)]], float)
    jacobian_inv = np.linalg.inv(jacobian)
    
    temp = jacobian_inv.dot(system_eq)
    
    temp = np.array([[x0_1], [x0_2]], float) - temp
    x0_1 = temp[0][0]
    x0_2 = temp[1][0]
    if abs(system_eq[0][0]) <= eps >= abs(system_eq[1][0]):
        print(x0_1, x0_2, i)
        break
        

eps = 10**-7
x0_1=1
x0_2=1

while True:
    system_eq = np.array([[4*(x0_1**2)-(x0_2**3)+28],[3*(x0_1**3)+4*(x0_2**2)-145]], float)
    
    jacobian = np.array([[8*x0_1, -3*(x0_2**2)],[9*(x0_1**2), 8*x0_2]], float)
    jacobian_inv = np.linalg.inv(jacobian)
    
    temp = jacobian_inv.dot(system_eq)
    
    temp = np.array([[x0_1], [x0_2]], float) - temp
    x0_1 = temp[0][0]
    x0_2 = temp[1][0]
    if abs(system_eq[0][0]) <= eps >= abs(system_eq[1][0]):
        print(x0_1, x0_2)
        break
        


def newton_nsae(x0_1, x0_2, eps, i=0):
    while True:
        i += 1
        system_eq = np.array([[2*(x0_1**3)+(x0_2**2)-3],[(10**(-1*x0_1))-x0_2-1.1]], float)

        jacobian = np.array([[6*(x0_1**2), 2*x0_2],[(-10**(-1*x0_1))*np.log(10), -1]], float)
        jacobian_inv = np.linalg.inv(jacobian)

        temp = jacobian_inv.dot(system_eq)

        temp = np.array([[x0_1], [x0_2]], float) - temp
        x0_1 = temp[0][0]
        x0_2 = temp[1][0]
        if abs(system_eq[0][0]) <= eps >= abs(system_eq[1][0]):
            return x0_1, x0_2, i, eps
            break
        
x0_1, x0_2 = float(input('Начальное приближение x = ')), float(input('Начальное приближение y = '))
x, y, iterations, eps = newton_nsae(x0_1, x0_2, 10 ** int(input('Степень точности = ')))

print()
print('Система уравнений: 2x^3+y^2=3')
print('\t\t   10^(-x)-y=1.1')
print('Точность:', eps)
print('Метод Ньютона')
print('x =', x, '|', 'y =', y)
print('n =', iterations)
print('Проверка:', 2*(x**3)+(y**2)-3, '|', (10**(-1*x))-y-1.1)


# Пример из методички
import matplotlib.pyplot as plt

c = np.linspace(-3, 4, 100)

y = np.sin(c-0.6)-1.6
x = (np.cos(c)/3)+0.3

fig, ax = plt.subplots()

ax.plot(c, y)
ax.set_xlim(-1,1)
ax.plot(x, c)
ax.grid()

fig.set_figwidth(15)
fig.set_figheight(8)

plt.show()

# Пример из методички

c = np.linspace(-3, 4, 100)

y = np.sin(c-0.6)-1.6
x = (np.cos(c)/3)+0.3

fig, ax = plt.subplots()

ax.plot(c, y)
ax.plot(x, c)
ax.grid()

fig.set_figwidth(15)
fig.set_figheight(8)

plt.show()


c = np.linspace(-5, 5, 150)

x = np.cbrt(4*(-1*((c)**2)+3))/2
y = (10**(-1*c))-1.1

fig, ax = plt.subplots()

ax.plot(x, c)
ax.plot(c, y)
ax.set_xlim(-2,4)
ax.set_ylim(-4,8)
ax.grid()

fig.set_figwidth(8)
fig.set_figheight(10)

plt.show()


c = np.linspace(-3, 3, 150)

x = np.cbrt(4*(-1*((c)**2)+3))/2
y = (10**(-1*c))-1.1

fig, ax = plt.subplots()

ax.plot(x, c)
ax.plot(c, y)
ax.set_xlim(-1,2)
ax.set_ylim(-1.5,2.5)
ax.grid()

fig.set_figwidth(8)
fig.set_figheight(10)

plt.show()
print('Из графика видно, что система имеет два решения.\nПервое решение в области D1: 0.8 < x < 1.2; -1.2 < y < -0.8')



#x_max = -0.1
#y_max = 2
x_max = 1.2
y_max = -0.8

d_f1_x = 0
d_f1_y = -1*(((2**(5/3))*y_max)/(6*((-1*(y_max**2)+3)**(2/3))))
d_f2_x = (-10**(-1*(x_max)))*np.log(10)
d_f2_y = 0

print(abs(d_f2_x), abs(d_f1_y))

if abs(d_f2_x) < 1 and abs(d_f1_y) < 1:
    print('|dф1/dx| + |dф2/dx| = ', abs(d_f2_x), '<', '1')
    print('|dф1/dy| + |dф2/dy| = ', abs(d_f1_y), '<', '1')
    print('Условия сходимости выполняются')
    
    
# Пример из методички
x0 = 0.15
y0 = -2
i = 0
while True:
    i += 1
    system_eq = np.array([[(np.cos(y0)/3)+0.3],[np.sin(x0-0.6)-1.6]], float)
    
    if abs(system_eq[0][0] - x0) and abs(system_eq[1][0] - y0) <= 10 ** -4:
        print(x0, y0, i)
        break
    x0 = system_eq[0][0]
    y0 = system_eq[1][0]
    
    

def fixed_iter(x0, y0, eps, i=0):
    while True:
        i += 1
        system_eq = np.array([[np.cbrt(4*(-1*((y0)**2)+3))/2],[(10**(-1*x0))-1.1]], float)

        if abs(system_eq[0][0] - x0) <= eps >= abs(system_eq[1][0] - y0):
            return x0, y0, i, eps
        x0 = system_eq[0][0]
        y0 = system_eq[1][0]
        
x0, y0 = float(input('Начальное приближение x = ')), float(input('Начальное приближение y = '))
x, y, iterations, eps = fixed_iter(x0, y0, 10 ** int(input('Степень точности = ')))

print()
print('Система уравнений: 2x^3+y^2=3')
print('\t\t   10^(-x)-y=1.1')
print('Точность:', eps)
print('Метод простой итерации')
print('x =', x, '|', 'y =', y)
print('n =', iterations)
print('Проверка:', 2*(x**3)+(y**2)-3, '|', (10**(-1*x))-y-1.1)
