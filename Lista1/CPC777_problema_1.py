# -*- coding: utf-8 -*-
"""
Created on Sat Sep 18 18:19:46 2021

@author: User
"""

#Resolução do Problema 1 - Lista 1:Fluxo radial unidimensional - analítico

import math as m
import matplotlib.pyplot as plt
import scipy.special as sc

#Dados de entrada
re = 60690                     #raio externo [cm]
rw = 8.9                       #raio do poço [cm]
h = 304.8                      #espessura do reservatório [cm]
pi = 149.7                     #pressão inicial [atm]
q0 = 5180                      #vazão [cm^3 std/ s]
k = 0.15                       #permeabilidade [Darcy]
phi= 0.2                       #porosidade [x100%]
ct = 0.0002204                 #compressibilidade total [1/atm]
mu = 0.33                      #viscosidade [cp]
Bo = 1.5                       # fato volume formação [m^3/m^3 std]
g = 1.78108                    #gamma para p em rw


C = Bo*q0*mu/(4*m.pi*k*h)      #parte constante da equação 3.225 de Rosa et al.

t1 = 60                         #[s] t = 1 min
t2 = 60*60                      #[s] t = 1 hora
t3 = 60*60*24*30                #[s] t = 30 dias 
t4 = 60*60*24*365               #[s] t = 1 ano
t5 = 60*60*24*365*10            #[s] t = 10 anos

#Definindo a lista de raios - r[0] = 0 metros. A correção será aplicada no cálculo da pressão.
r = []                                                     #lista a ser montada
n = 1000                                                   #quantidade de termos  

for i in range(n+1):
    r = r + [i*(re-rw)/n]
    
print(r[-1])                                               #teste para saber o último valor: r[n+1]

#Prerfil de pressão para t = 1min
p1 = []                                                    #lista a ser montada                                         

for i in range(n+1):
    p1 = p1 + [pi+C*sc.expi(-phi*mu*ct*r[i]**2/(4*k*t1))]  #montagem da lista para p1@1min

p1[0]= pi+C*sc.expi(-phi*mu*ct*rw**2/(4*k*t1))             #correção para r[0] = rw  

print(p1[-1])                                              #teste para saber o último valor: p1[n+1]

   
#Prerfil de pressão para t = 1hora
p2 = []

for i in range(n+1):
    p2 = p2 + [pi+C*sc.expi(-phi*mu*ct*r[i]**2/(4*k*t2))]

p2[0]= pi+C*sc.expi(-phi*mu*ct*rw**2/(4*k*t2)) 

print(p2[-1])


#Prerfil de pressão para t = 30 dias
p3 = []

for i in range(n+1):
    p3 = p3 + [pi+C*sc.expi(-phi*mu*ct*r[i]**2/(4*k*t3))]

p3[0]= pi+C*sc.expi(-phi*mu*ct*rw**2/(4*k*t3)) 

print(p3[-1])

#Prerfil de pressão para t = 1 ano
p4 = []

for i in range(n+1):
    p4 = p4 + [pi+C*sc.expi(-phi*mu*ct*r[i]**2/(4*k*t4))]

p4[0]= pi+C*sc.expi(-phi*mu*ct*rw**2/(4*k*t4)) 

print(p4[-1])

#Prerfil de pressão para t = 10 anos
p5 = []

for i in range(n+1):
    p5 = p5 + [pi+C*sc.expi(-phi*mu*ct*r[i]**2/(4*k*t5))]

p5[0]= pi+C*sc.expi(-phi*mu*ct*rw**2/(4*k*t5)) 

print(p5[-1])

#Gráficos
plt.plot(r,p1, label="1 minuto")
plt.plot(r,p2, label="1 hora")
plt.plot(r,p3, label="30 dias")
plt.plot(r,p4, label="1 ano")
plt.plot(r,p5, label="10 anos")
plt.xlabel('raio [m]')
plt.ylabel('Pressão [atm]')
plt.legend()
plt.show()
