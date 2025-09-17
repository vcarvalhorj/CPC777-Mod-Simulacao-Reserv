# -*- coding: utf-8 -*-
"""
Created on Sun Sep 26 10:40:38 2021

@author: vcarv
"""

#Resolução do Problema 2 - Lista 2:Fluxo radial unidimensional - numérico
#Método explícito

import math as m
import scipy.special as sc
import matplotlib.pyplot as plt


#1 - Dados de entrada
re = 60690                                                             #raio externo [cm]
rw = 8.9                                                               #raio do poço [cm]
h = 304.8                                                              #espessura do reservatório [cm]
pi = 149.7                                                             #pressão inicial [atm]
q0 = 5180                                                              #vazão [cm^3 std/ s]
k = 0.15                                                               #permeabilidade [Darcy]
phi= 0.2                                                               #porosidade [x100%]
ct = 0.0002204                                                         #compressibilidade total [1/atm]
mu = 0.33                                                              #viscosidade [cp]
Bo = 1.5                                                               # fato volume formação [m^3/m^3 std]
g = 1.78108                                                            #gamma para p em rw


#2 - Cálculo das constantes

dt = 2                                                                 #passo de tempo [s]
dr = 303.45                                                            #espassamento da malha [cm]

ni = k/(phi*mu*ct)                                                     #constante de difusividade hidráulica - eq.3.113 de Rosa et al.

#parcelas constantes da transmissibilidade do nó
tau_c1 = ni*dt/(dr**2)          
tau_c2 = ni*dt/(2*dr)                                                  #parcela que multiplica 1/ri no cálculo da pressão

#parcelas constantes B - relativo às pressões no poço
cc = Bo*q0*mu/(2*m.pi*k*h)                                             


#3 - Definição da lista de raios - pontos i - que serão calculados o campos de pressão

r = []
n = int (re/dr)
for i in range(n+1):
    r += [i*dr]
r[0] = rw                                                            #correção r1 = 0 para r1 = rw

print(r[-1])                                                         #teste para saber o último valor: r[n+1]

Nr = len(r)                                                          #numero de pontos no espaço


#4 - Condições iniciais

t1 = 60*60*24*365                                                    #[s] tempo total de análise
Nt = int(t1/dt)                                                         #número de passos no tempo

Pi = []                                                               #campo de pressões iniciais
for i in range(Nr):                                                   #Pi = [pi for i in range(Nr)]
    Pi += [pi]

    
#5 - Cálculo do campos de pressões no tempo (numérico)

P = []                                                                #pressões no tempo presente
P_aux = []                                                            #pressões no tempo futuro
for i in range (Nt):
    if i == 0:
        P = Pi.copy()
        print('it =', i,'// t[s] =', round(i*dt, 2), '// horas =', round(i*dt/3600, 2), ' // dias =', round(i*dt/3600/24, 2), 
              '// mes = ', round(i*dt/3600/24/30.2, 1), '// ano =', round(i*dt/3600/24/365.5, 1), '// pw[atm] =',round(P[0], 2))
    else:
        P = P_aux.copy()
        P_aux = []
        print('it =', i,'// t[s] =', round(i*dt, 2), '// horas =', round(i*dt/3600, 2), ' // dias =', round(i*dt/3600/24, 2), 
              '// mes = ', round(i*dt/3600/24/30.2, 1), '// ano =', round(i*dt/3600/24/365.5, 1), '// pw[atm] =',round(P[0], 2))
    for j in range (Nr):
        if j == 0:
            P_aux += [4*tau_c1*P[j+1]+(1-4*tau_c1)*P[j]-cc]                
        elif j == Nr-1:
            P_aux += [P[j]-cc/(r[j]*dr)]
        else:
            P_aux += [((1/r[j])*tau_c2 +tau_c1)*P[j+1] + (tau_c1-(1/r[j])*tau_c2)*P[j-1] + (1 - 2*tau_c1)*P[j] -cc/r[j]]
    

print (P_aux[-1])



#6 - Perfil de pressão  (analítico)
C = Bo*q0*mu/(4*m.pi*k*h)                                           #parte constante da equação 3.225 de Rosa et al.

p1 = []                                                             #lista a ser montada                                         

for i in range(Nr):
    p1 = p1 + [pi+C*sc.expi(-phi*mu*ct*r[i]**2/(4*k*t1))]           #montagem da lista para p1@t1


print(p1[-1])                                                       #teste para saber o último valor: p1[n+1]


#7 - Gráficos
plt.plot(r, P_aux, label="Numérico")
plt.plot(r, p1, label="Analítico")
plt.xlabel('raio [cm]')
plt.ylabel('Pressão [atm] - 1 ano')
plt.legend()
plt.show()






