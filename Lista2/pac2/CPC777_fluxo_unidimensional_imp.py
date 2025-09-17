# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 15:50:37 2021

@author: User
"""

#Resolução do Problema 2 - Lista 2:Fluxo radial unidimensional - numérico
#Método implícito

import math as m
import scipy.special as sc
import matplotlib.pyplot as plt

#1 - Dados de entrada
re = 60690                                                             #raio externo [cm]
rw = 8.9                                                               #raio do poço [cm]
h = 304.8                                                              #espessura do reservatório [cm]
pi = 149.7                                                             #pressão inicial [atm]
q0 = 5180                                                              #vazão [cm^3 std/ s]
kk = 0.15                                                              #permeabilidade [Darcy]
phi= 0.2                                                               #porosidade [x100%]
ct = 0.0002204                                                         #compressibilidade total [1/atm]
mu = 0.33                                                              #viscosidade [cp]
Bo = 1.5                                                               # fato volume formação [m^3/m^3 std]
g = 1.78108                                                            #gamma para p em rw


#2 - Cálculo das constantes

dt = 5                                                              #passo de tempo [s]
dr = 121.38                                                           #espassamento da malha [cm]

ni = kk/(phi*mu*ct)                                                   #constante de difusividade hidráulica - eq.3.113 de Rosa et al.

#parcelas constantes da transmissibilidade do nó
tau_c1 = ni*dt/(dr**2)          
tau_c2 = ni*dt/(2*dr)                                                 #parcela que multiplica 1/ri no cálculo da pressão


# parcelas constantes relativo ao pw 
cc = Bo*q0*mu/(2*m.pi*kk*h)    

#3 - Definição da lista de raios - pontos i - que serão calculados o campos de pressão
#A quantidade de nós, Nr, irá definir o tamanho da matriz das transmissibiliades - Nr x Nr

r = []
n = int (re/dr)
for i in range(n+1):
    r += [i*dr]
r[0] = rw                                                            #correção r1 = 0 para r1 = rw

print(r[-1])                                                         #teste para saber o último valor: r[n+1]

Nr = len(r)                                                          #numero de pontos no espaço


#4 - Condições iniciais

t1 = 60*60*24                                                         #[s] tempo total de análise
Nt = int(t1/dt)                                                       #número de passos no tempo

Pi = []                                                              #campo de pressões iniciais
for i in range(Nr):                                                  #Pi = [pi for i in range(Nr)]
    Pi += [pi]
    
   
#5 - Montagem da matriz de transmissibilidades

#matriz zerada
tr = []
for i in range (Nr):
    tr += [[]]
for i in range (Nr):
    tr[i] += [0 for j in range(Nr)]                                  #matriz = [[0 for i in range(Nr)] for j in range (Nr)].Montagem de matriz - lista de lista

#matriz tri - diagonal                                 
for i in range (Nr):
    for j in range (Nr):
        if i==j:                                                      #parcela diagonal
            tr [i][j] = 1 + 2*tau_c1                              
        elif j-i == 1:                                                #diagonal superior
            if i == 0:
                tr [i][j] = -2*tau_c1
            else:
                tr [i][j] = -(tau_c2*(1/r[i]) + tau_c1)                
        elif j-i == -1:                                               #diagonal inferior
             if i == (Nr-1):
                tr [i][j] = -2*tau_c1
             else:
                tr [i][j] = -(tau_c1 - tau_c2*(1/r[i]))                



#6 - Inversão da matriz de transmissibiliade                    
#https://numpy.org/doc/stable/reference/generated/numpy.linalg.inv.html?highlight=inv#numpy.linalg.inv


from numpy.linalg import inv
tr_inv = inv(tr)
      
#7 - Cálculo do campos de pressões no tempo (numérico)   
    
P = []                                                                           #pressões no tempo presente
P_aux = []                                                                       #pressões no tempo futuro

for k in range (Nt):
    if k == 0:
        P = Pi.copy()
        print('it =', k,'// t[s] =', round(k*dt, 2), '// horas =', round(k*dt/3600, 2), ' // dias =', round(k*dt/3600/24, 2), 
              '// mes = ', round(k*dt/3600/24/30.2, 1), '// ano =', round(k*dt/3600/24/365.5, 1), '// pw[atm] =',round(P[0], 2))
    else:
        P = P_aux.copy()
        P_aux = []
        print('it =', k,'// t[s] =', round(k*dt, 2), '// horas =', round(k*dt/3600, 2), ' // dias =', round(k*dt/3600/24, 2), 
              '// mes = ', round(k*dt/3600/24/30.2, 1), '// ano =', round(k*dt/3600/24/365.5, 1), '// pw[atm] =',round(P[0], 2))
    for i in range (Nr):
        P_aux += [0]
        if i == 0:
            for j in range (Nr):
                if j == 0:
                    P_aux [i] = P_aux[i] + tr_inv[i][j]*(P[j] - (ni*dt/(r[i]**2))*cc)
                else:
                    P_aux [i] = P_aux[i] + tr_inv[i][j]*P[j] 
        else:
            for j in range (Nr):
                P_aux [i] = P_aux[i] + tr_inv[i][j]*P[j] 
            
print (P_aux[-1])     
        
#8 - Perfil de pressão (analítico)
C = Bo*q0*mu/(4*m.pi*kk*h)                                           #parte constante da equação 3.225 de Rosa et al.

p1 = []                                                              #lista a ser montada                                         

for i in range(Nr):
    p1 = p1 + [pi+C*sc.expi(-phi*mu*ct*r[i]**2/(4*kk*t1))]           #montagem da lista para p1@t1


print(p1[-1])                                                        #teste para saber o último valor: p1[n+1]


#7 - Gráficos
plt.plot(r, P_aux, label="Numérico")
plt.plot(r, p1, label="Analítico")
plt.xlabel('raio [cm]')
plt.ylabel('Pressão [atm] - 1 dia')
plt.legend()
plt.show()
        
        
        


    
    
    
    