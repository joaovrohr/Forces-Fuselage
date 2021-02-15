import numpy as np
import matplotlib.pyplot as plt

Cn_a_5_aletas = 4.74      #[1/rad] Coeficiente de força normal derivado com relação ao ângulo de ataque para as aletas
Cn_a_5_motor = 0.614      #[1/rad] Coeficiente de força normal derivado com relação ao ângulo de ataque para o motor
Cn_a_5_tanque = 0.653     #[1/rad] Coeficiente de força normal derivado com relação ao ângulo de ataque para o tanque
Cn_a_5_avionica = 0.312   #[1/rad] Coeficiente de força normal derivado com relação ao ângulo de ataque para a aviônica
Cn_a_5_rec = 0.429        #[1/rad] Coeficiente de força normal derivado com relação ao ângulo de ataque para a recuperação
Cn_a_5_coifa = 2.16       #[1/rad] Coeficiente de força normal derivado com relação ao ângulo de ataque para a coifa

Cn_a_5 = [Cn_a_5_aletas, Cn_a_5_motor, Cn_a_5_tanque, Cn_a_5_avionica, Cn_a_5_rec, Cn_a_5_coifa]

m_aletas = 0.101          #[kg] Massa da seção de aletas
m_motor = 6.0933          #[kg] Massa da seção de motor
m_tanque = 6.418          #[kg] Massa da seção de tanque
m_avionica = 3.8015       #[kg] Massa da seção de aviônica
m_rec = 1.6665            #[kg] Massa da seção de recuperação
m_coifa = 0.5155          #[kg] Massa da seção de coifa

m = [m_aletas, m_motor, m_tanque, m_avionica, m_rec, m_coifa]

M_total = m_aletas + m_rec + m_avionica + m_tanque + m_motor + m_coifa  #[kg] Massa total
alfa = 5*np.pi/180                                                      #[rad] Ângulo de ataque
beta = 5*np.pi/180                                                      #[rad] Ângulo de direção de voo

S = 0.0123                   #[m^2] Área de referência
vel = 150                      #[m/s] Velocidade máxima
rho_ar = 1.225               #[kg/m^3] Densidade do ar
q = 0.5*rho_ar*(vel**2)        #[Pa] Pressão dinâmica

D_e = 0.125                            #[m] Diâmetro externo da fuselagem
D_i = 0.121                            #[m] Diâmetro interno da fuselagem
At = ((D_e/2)**2 - (D_i/2)**2)*np.pi   #[m^2] Área da seção transversal da fuselagem
rho_fibra = 1850                       #[kg/m^3] Densidade da fibra de vidro

CG = 1.09                              #[m] CG na pressão dinâmica máxima
L = 2.06                               #[m] Comprimento da fuselagem

#Calculando as forças normais aerodinâmicas

CN = 0         #[N] Força aerodinâmica normal total
Cn = []        #[N] Vetor com as forças aerodinâmicas normais para cada seção

for i in range(len(Cn_a_5)):
    Cn.append(q*S*alfa*Cn_a_5[i])
    CN += q*S*alfa*Cn_a_5[i]

#Calculando as forças peso

P = 0         #[N] Força peso total
P_vec = []    #[N] Vetor com as forças peso

for i in range(len(m)):
    P_vec.append(9.81*m[i]*np.sin(alfa+beta))
    P += 9.81*m[i]*np.sin(alfa+beta)

#Calculando a aceleração normal

ay = (CN-P)/M_total       #[m/s^2] Aceleração normal

#Locais de aplicação de forças
CP_aletas = 1.98                    #[m] CP da seção das aletas (0 na junção da fuselagem com a coifa)
CP_motor = 1.75                     #[m] CP da seção do motor (0 na junção da fuselagem com a coifa)
CP_tanque = 1.09                    #[m] CP da seção do tanque (0 na junção da fuselagem com a coifa)
CP_avionica = 0.6                   #[m] CP da seção da aviônica (0 na junção da fuselagem com a coifa)
CP_rec = 0.22                       #[m] CP da seção da recuperação (0 na junção da fuselagem com a coifa)
CP_coifa = 0                        #[m] CP da seção da coifa (0 na junção da fuselagem com a coifa)

CG_coifa = 0.154                    #[m] CG da seção de coifa (utilizado para o cálculo da força inercial angular da coifa)

CP_vec = [CP_aletas, CP_motor, CP_tanque, CP_avionica, CP_rec, CP_coifa]

x1 = 2.06                           #[m] Posição do primeiro bulkhead o motor (perto da tubeira) (0 na junção da fuselagem com a coifa)
x2 = 2.00                           #[m] Posição do CG das aletas (0 na junção da fuselagem com a coifa)
x3 = 1.42                           #[m] Posição da junção do motor com o tanque (0 na junção da fuselagem com a coifa)
x4 = 0.75                           #[m] Posição da junção do tanque com a aviônica (0 na junção da fuselagem com a coifa)
x5 = 0.43                           #[m] Posição da junção da aviônica com a recuperação (0 na junção da fuselagem com a coifa)
x6 = 0                              #[m] Posição da junção da junção da recuperação com a coifa (0 na junção da fuselagem com a coifa)

x_vec = []              #[m] Vetor que guarda os valores apresentados logo acima

#Discretização da fuselagem

n_inter = 1000          #[-] Número de pontos na discretização
step = L/n_inter        #[m] Tamanho do intervalo entre cada ponto
x = []                  #[m] Vetor de discretização dos pontos ao longo da fuselagem

for i in range(n_inter+1):
    x.append(i*step)

#Calculando aceleração angular

I = 10.6                #[kg.m^2] Momento de inércia de massa do foguete
M = 0                   #[Nm] Momento total agindo no foguete

for i in range(len(Cn)):
    M += (CG-CP_vec[i])*Cn[i]

R = M/I                 #[rad/s^2] Aceleração angular resultante do foguete

#Calculando os esforços cortantes

Vs = []         #[N] Vetor que guarda os esforços cortantes
e1 = 0          #variável ativadora de seção
e2 = 0          #variável ativadora de seção
e3 = 0          #variável ativadora de seção
e4 = 0          #variável ativadora de seção
e5 = 0          #variável ativadora de seção
e6 = 0          #variável ativadora de seção
e7 = 0          #variável ativadora de seção
e8 = 0          #variável ativadora de seção
e9 = 0          #variável ativadora de seção
e10 = 0         #variável ativadora de seção
e11 = 0         #variável ativadora de seção

#Aqui se calcula os esforços cortantes, a ideia é prosseguir em x (Como fazer o gráfico dos esforços cortantes) e se x
#passar por uma força ela será levada em consideração pelas variáveis ativadoras.

#Começa-se pela parte traseira do foguete.

for i in range(len(x)):
    if 0 < x[i] and x[i] < L-x2:
        e1 = 1
    if L-x2 <= x[i] and x[i] < L-CP_aletas:
        e2 = 1
    elif L-CP_aletas <= x[i] and x[i] < L-CP_motor:
        e3 = 1
    elif L-CP_motor <= x[i] and x[i] < L-x3 :
        e4 = 1
    elif L-x3 <= x[i] and x[i] < L-CP_tanque:
        e5 = 1
    elif L-CP_tanque <= x[i] and x[i] < L-x4:
        e6 = 1
    elif L-x4 <= x[i] and x[i] < L-CP_avionica:
        e7 = 1
    elif L-CP_avionica <= x[i] and x[i] < L-x5:
        e8 = 1
    elif L-x5 <= x[i] and x[i] < L-CP_rec:
        e9 = 1
    elif L-CP_rec <= x[i] and x[i] < L:
        e10 = 1
    elif x[i] == L:
        e11 = 1
    Vs.append(-(P_vec[1]/2)*e1 - (P_vec[0])*e2 + (Cn[0])*e3 + (Cn[1])*e4 - (P_vec[1]/2 + P_vec[2]/2)*e5 + (Cn[2])*e6 - (P_vec[2]/2+P_vec[3]/2)*e7 +
              (Cn[3])*e8 - (P_vec[3]/2 + P_vec[4]/2)*e9 + (Cn[4])*e10 + (Cn[5] - P_vec[5] - P_vec[4]/2)*e11 + R*(rho_fibra*At*x[i]*(x[i]/2-CG) + (m_motor/2)*e1*(x1-CG) +
              (m_aletas)*e2*(x2-CG) + (m_motor/2 + m_tanque/2)*e5*(x3-CG) + (m_tanque/2+m_avionica/2)*e7*(x4-CG) +
              (m_avionica/2 + m_rec/2)*e9*(x5-CG) + ((m_rec/2)*(x6-CG) + m_coifa*(CG_coifa-(CG+0.25)))*e11) -
              ay*(rho_fibra*At*x[i] + (m_motor/2)*e1 + (m_aletas)*e2 + (m_motor/2 + m_tanque/2)*e5 + (m_tanque/2+m_avionica/2)*e7 +
              (m_avionica/2 + m_rec/2)*e9 + (m_rec/2 + m_coifa)*e11))


#Calculando o momento fletor

Mb = []         #[Nm] Vetor que guarda os momentos fletores
MB = 0          #[Nm] Momento fletor acumulado

#Aqui se faz uma integração pela regra do trapézio
for i in range(len(x)-1):
    Mb.append(MB + (Vs[i+1]+Vs[i])*step/2)
    MB += (Vs[i+1]+Vs[i])*step/2

Vs[len(Vs)-1] = 0 #Ajeitada no gráfico
Vs[0] = 0         #Ajeitada também
Mb.append(0)      #Outra ajeitada

Vsmax = 0         #[N] Esforço cortante máximo
pos1 = 0          #[m] Posição do esforço cortante máximo

for i in range(len(Vs)): #Busca pelo esforço cortante máximo
    if abs(Vs[i]) > Vsmax:
        Vsmax = Vs[i]
        pos1 = x[i]

print("Esforço cortante máximo:", Vsmax, "Em x =", pos1)

Mbmax = 0        #[Nm] Momento fletor máximo
pos2 = 0         #[m] Posição do momento fletor máximo

for i in range(len(Mb)): #Busca pelo momento fletor máximo
    if abs(Mb[i]) > Mbmax:
        Mbmax = Mb[i]
        pos2 = x[i]

print("Momento fletor máximo:", Mbmax, "Em x =", pos2)

#Plot da força cortante
plt.plot(x, Vs, "-")                    #Passando os pontos para plotar
plt.xlabel("Posição (m)")               #Dando nome ao eixo x
plt.ylabel("Esforço Cortante V (N)")    #Dando nome ao eixo y
plt.plot(pos1, Vsmax, ".", color = 'r', label = "$V_{máx}$ = " + "{:.3f}".format(Vsmax) + " N")  #Passando a posição e nota de ponto máximo
plt.legend()                            #Comando que mostra as legendas
plt.show()                              #Comando que mostra o gráfico

#Plot do momento fletor
plt.plot(x, Mb, "-")                    #Passando os pontos para plotar
plt.xlabel("Posição (m)")               #Dando nome ao eixo x
plt.ylabel("Momento Fletor M (N.m)")    #Dando nome ao eixo y
plt.plot(pos2, Mbmax, ".", color = 'r', label = "$M_{máx}$ = " + "{:.3f}".format(Mbmax) + " Nm") #Passando a posição e nota de ponto máximo
plt.legend()                            #Comando que mostra as legendas
plt.show()                              #Comando que mostra o gráfico

E = 7*(10**9)                           #[Pa] Módulo de elasticidade da fibra de vidro
Iz = np.pi*((D_e**4) - (D_i**4))/64     #[m^4] Momento de inércia de área da seção transversal da fuselagem

dvdx = []                               #Rotação da fuselagem (eu acho)
DVDX = 0                                 #Rotação acumulada (eu acho também)

for i in range(len(Mb)-1):               #Integrando pelo método do trapézio também
    dvdx.append(DVDX + (Mb[i+1]+Mb[i])*step/2)
    DVDX += (Mb[i+1]+Mb[i])*step/2

dvdx.append(dvdx[len(dvdx)-1])           #Adicionando ponto no gráfico de novo

v = []                     #[mm] Vetor com os deslocamentos verticais da fuselagem
V = 0                      #[mm] Deslocamento acumulado da fuselagem

for i in range(len(dvdx)-1):             #Integração pela regra dos trapézios
    v.append(V + (dvdx[i+1]+dvdx[i])*step*1000/(2*E*Iz))
    V += (dvdx[i+1]+dvdx[i])*step*1000/(2*E*Iz)

v.append(v[len(v)-1])     #Ajeitando os pontos

print("Deformação vertical máxima:", v[len(v)-1])

#Plotando as deformações
plt.plot(x, v, "-")                        #Passando os pontos para plotar
plt.xlabel("Posição (m)")                  #Dando nome ao eixo x
plt.ylabel("Deformação Vertical v (mm)")   #Dando nome ao eixo y
plt.plot(x[len(x)-1], v[len(v)-1], ".", color = 'r', label = "$v_{máx}$ = " + "{:.3f}".format(v[len(v)-1]) + " mm")  #Passando a posição e nota de ponto máximo
plt.legend()                               #Comando que mostra as legendas
plt.show()                                 #Comando que mostra o gráfico

Q = ((D_e**3)-(D_i)**3)/12                 #[m^3] Momento estático de inércia
t = D_e-D_i                                #[m] Duas vezes a espessura da parede

tau_max = Vsmax*Q/(Iz*t)                   #[Pa] Tensão cisalhante na seção transversal da fuselagem

print("Tensão de cisalhamento máxima:", tau_max)

sigma = Mbmax*D_e/(Iz*2)                   #[Pa] Tensão normal máxima gerada pelo momento fletor na seção

print("Tensão normal máxima gerada pelo momento fletor:", sigma)

#Cálculo das forças axiais

Cd_aletas = 0.05          #Coeficiente de arrasto das aletas
Cd_motor = 0.23           #Coeficiente de arrasto do motor
Cd_tanque = 0.08          #Coeficiente de arrasto do tanque
Cd_avionica = 0.04        #Coeficiente de arrasto da aviônica
Cd_recuperacao = 0.05     #Coeficiente de arrasto da recuperação
Cd_coifa = 0.02           #Coeficiente de arrasto da coifa

Cd_vec = [Cd_aletas, Cd_motor, Cd_tanque, Cd_avionica, Cd_recuperacao, Cd_coifa]

Fd_aletas = Cd_aletas*S*q              #Força de arrasto no CP das aletas
Fd_motor = Cd_motor*S*q                #Força de arrasto no CP do motor
Fd_tanque = Cd_tanque*S*q              #Força de arrasto no CP do tanque
Fd_avionica = Cd_avionica*S*q          #Força de arrasto no CP da aviônica
Fd_recuperacao = Cd_recuperacao*S*q    #Força de arrasto no CP da recuperação
Fd_coifa = Cd_coifa*S*q                #Força de arrasto no CP da coifa

Fd_vec = [Fd_aletas, Fd_motor, Fd_tanque, Fd_avionica, Fd_recuperacao, Fd_coifa]

T = 1800        #[N] Empuxo

#Cálculo da aceleração axial

FD = 0
for i in range(len(Fd_vec)-1):
    FD += Fd_vec[i]
ax = (T - FD - 9.81*M_total)/M_total

D = []          #Vetor que guarda os esforços axiais
e1 = 0          #variável ativadora de seção
e2 = 0          #variável ativadora de seção
e3 = 0          #variável ativadora de seção
e4 = 0          #variável ativadora de seção
e5 = 0          #variável ativadora de seção
e6 = 0          #variável ativadora de seção
e7 = 0          #variável ativadora de seção
e8 = 0          #variável ativadora de seção
e9 = 0          #variável ativadora de seção
e10 = 0         #variável ativadora de seção
e11 = 0         #variável ativadora de seção

for i in range(len(x)):
    if 0 < x[i] and x[i] < L-x2:
        e1 = 1
    if L-x2 <= x[i] and x[i] < L-CP_aletas:
        e2 = 1
    elif L-CP_aletas <= x[i] and x[i] < L-CP_motor:
        e3 = 1
    elif L-CP_motor <= x[i] and x[i] < L-x3 :
        e4 = 1
    elif L-x3 <= x[i] and x[i] < L-CP_tanque:
        e5 = 1
    elif L-CP_tanque <= x[i] and x[i] < L-x4:
        e6 = 1
    elif L-x4 <= x[i] and x[i] < L-CP_avionica:
        e7 = 1
    elif L-CP_avionica <= x[i] and x[i] < L-x5:
        e8 = 1
    elif L-x5 <= x[i] and x[i] < L-CP_rec:
        e9 = 1
    elif L-CP_rec <= x[i] and x[i] < L:
        e10 = 1
    elif x[i] == L:
        e11 = 1
    D.append(T*e1/2 - (9.81 + ax)*(m_motor/2)*e1 - (9.81 + ax)*(m_aletas/2)*e2 - Fd_aletas*e3 - Fd_motor*e4 +
             (9.81 + ax)*(m_motor/2 + m_tanque/2)*e5 + T*e5/2 - Fd_tanque*e6 - (9.81 + ax)*(m_avionica/2 + m_tanque/2)*e7 -
             Fd_avionica*e8 - (9.81 + ax)*(m_avionica/2 + m_rec/2)*e9 - Fd_recuperacao*e10 - (9.81 + ax)*(m_rec/2 + m_coifa)*e11 -
             Fd_coifa*e11 - rho_fibra*At*x[i]*(ax + 9.81))

D_max = 0                                     #Força axial no ponto de máximo momento fletor (bulkhead do tanque com o motor)
pos3 = 0                                      #Posição da força axial máxima

for i in range(len(D)):
    if abs(D[i]) > D_max:
        D_max = abs(D[i])
        pos3 = x[i]

print("Esforço axial máximo:", D_max, "Em x =", pos3)

#Plot da força axial
plt.plot(x, D, "-")                    #Passando os pontos para plotar
plt.xlabel("Posição (m)")               #Dando nome ao eixo x
plt.ylabel("Esforço Axial D (N)")    #Dando nome ao eixo y
plt.plot(pos3, D_max, ".", color = 'r', label = "$D_{máx}$ = " + "{:.3f}".format(D_max) + " N")  #Passando a posição e nota de ponto máximo
plt.legend()                            #Comando que mostra as legendas
plt.show()                              #Comando que mostra o gráfico

sigma2 = D_max*4/(((D_e**2)-(D_i**2))*np.pi)  #[Pa] Tensão normal máxima gerada pela força axial máxima

print("Tensão normal máxima gerada pela força axial:", sigma2)

print("Tensão normal máxima total:", sigma+sigma2)

z = np.linspace(-D_e/2, D_e/2, 300) #Criando vetor de pontos de posição vertical na fuselagem

plt.plot(z, (sigma2 - Mbmax*z/(Iz))/1000000, "-")    #Passando os pontos para plotar
plt.xlabel("Posição vertical na seção (m)")          #Dando nome ao eixo x
plt.ylabel("Tensão Normal (MPa)")                    #Dando nome ao eixo y
plt.plot(z[0], (sigma2 - Mbmax*z[0]/(Iz))/1000000, ".", color = 'r', label = "$\sigma_{máx}$ = " + "{:.3f}".format((sigma2 - Mbmax*z[0]/(Iz))/1000000) + " MPa")
plt.grid()                                           #Comando que mostra um grid
plt.legend()                                         #Comando que mostra as legendas
plt.show()                                           #Comando que mostra o gráfico