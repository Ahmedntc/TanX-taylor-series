import math
import numpy as np
import matplotlib.pyplot as plt




#diff aproach
def bernoulli(z):
    if z == 0:
        return 1
    else:
        total = 0
        for i in range (0, z + 1):
            #math.comb nos retorna o coeficiente binomial de uma expansao que no caso da tangente é o bernoulli
            total += math.comb(z, i) * bernoulli(i) / (z - i + 1)
        return 1 - total

def expansaoTanX(x0, n):
    somaTanx = 0
    for i in range(1, n+1):
        somaTanx += ((bernoulli(2 * i)) / math.factorial(2 * i)) * ((-4)**i) * (1 - (4**i)) * (x0**(2 * i - 1))
    return somaTanx


#aproximaçao da tangente por taylor / maclaurin utilizando numero de bernoulli
def tan_Taylor(x0, n):
    somaTan = 0
    for i in range(1 , n + 1):
    #calculando numero de bernoulli que sao os coeficientes da serie de taylor de tan x
        B = 0
        numB = 2 * i
        for j in range (numB + 1):
            aux = 0
            for k in range(0, j + 1):
                aux += pow(-1, k) * math.factorial(j)  * pow(k, numB) / (math.factorial(k) * math.factorial(j - k))   

            B += aux / ((j + 1))

        somaTan +=  pow(-4, i) * (1 - pow(4, i)) * B * pow(x0, 2 * i - 1) / math.factorial(2 * i)
    return somaTan

#aproximaçao da tangente usando relações trigonométricas
def tanSinCos(x0, n):
    tan = np.divide(seno(x0, n), cosseno(x0, n))
    return tan

#aproximação de cosseno usando taylor / maclaurin
def cosseno(x0, n):
  numerador = -np.multiply(x0,x0)
  denominador = 2.0
  cosseno = np.add(1.0, numerador/denominador) 

  for n in range(2, n):
    numerador = np.multiply(-numerador,np.multiply(x0,x0))
    denominador = 2*np.multiply(denominador,2*np.power(n,2)-n)
    cosseno = np.add(cosseno, np.divide(numerador,denominador))

  return cosseno

#aproximação de seno usando taylor / maclaurin
def seno(x0, n):
    seno = 0
    for i in range(n):
        sinal = (-1)**i  
        seno +=(x0**(2.0*i + 1))/math.factorial(2*i+1)*sinal
    return seno 



print("Tangente numpy: %.20f" % np.tan(1))
tanExpansion = tan_Taylor(1, 6)
print("Tangente aproximada por taylor/maclaurin: %.20f" %tanExpansion)
tanExpansion2 = expansaoTanX(1, 1)
print("Tangente aproximada por taylor/maclaurin versao 2: %.20f" %tanExpansion2)

tanSC = tanSinCos(1, 20)
print("Tangente aproximada por taylor/maclaurin seno cossseno: %.20f" %tanSC)





x = np.linspace(-2*np.pi, 2*np.pi, 100)
y1 = np.tan(1)
y2 = tanSinCos(1, 20)

fig, ax = plt.subplots()
fig.set_size_inches(20, 8)
ax.plot(1 , y1, linewidth = 2.0, label ="numpy Tan")
ax.plot(1 , y2, linewidth = 2.0, label ="Tan")

ax.legend()
ax.grid()

plt.savefig("fig.png")