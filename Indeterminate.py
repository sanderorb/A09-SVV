#### Imports ####

import numpy as np
from math import *
import matplotlib.pyplot as plt
from scipy.sparse import linalg

#### Definitions ####

def MCLY(F1,F2,x):
    My= F1*toggle(x,x1)**1
    return My
    
def toggle(x,loc):
    if x>=loc:
        return (x-loc)
    else:
        return 0
def  toggle0(x,loc):
    if x>=loc:
        return (1)
    else:
        return 0

def integrate(f,a,b):
    x = np.linspace(a,b,len(f)) # N+1 points make N subintervals
    y = f
    
    dx = x[1]-x[0]
    Int=[]
    for i in range (0,len(f)-1):
        T = (dx/2) * (y[i+1] + y[i])
        if i==0:
            Int.append(T)
        else:
             Int.append(T+Int[i-1])
    return Int

def Coefs_equi(x_1, x_2, x_3, x_a, s, h_a, l_a, theta, P, q_x, t_x):
    x_I = x_2 - (x_a / 2)
    x_II = x_2 + (x_a / 2)

    # x=[F_z1,F_z2,F_z3,F_y1,F_y2,F_y3,H,C1,C2,C3,C4,C5]
    M[0] = [ -(l_a - x_1) , -(l_a - x_2) , -(l_a - x_3) , 0 , 0 , 0 , cos(theta) * (l_a - x_I) , 0 , 0 , 0 , 0 , 0 ]
    M[1] = [ 0 , 0 , 0 , (l_a - x_1) , (l_a - x_2) , (l_a - x_3) , -sin(theta) * (l_a - x_I) , 0 , 0 , 0 , 0 , 0 ]
    M[2] = [ 0 , 0 , 0 , -(s - h_a/2) , -(s - h_a/2) , -(s - h_a/2) , -(h_a/2 * cos(theta) - s * sin(theta)) , 0 , 0 , 0 , 0 , 0 ]
    M[3] = [ 0 , 0 , 0 , 1 , 1 , 1 , -sin(theta) , 0 , 0 , 0 , 0 , 0 ]
    M[4] = [ -1 , -1 , -1 , 0 , 0 , 0 , cos(theta) , 0 , 0 , 0 , 0 , 0 ]

    b[0] = -P * cos(theta) * (l_a - x_II)
    b[1] = P * sin(theta) * (l_a - x_II) + integrate(integrate(q_x, 0, l_a), 0 , l_a)[-1]
    b[2] = P * (h_a/2 * cos(theta) - s * sin(theta)) - integrate(t_x, 0, l_a)[-1]
    b[3] = P * sin(theta) + integrate(q_x, 0, l_a)[-1]
    b[4] = - P * cos(theta)
    return M[0], M[1], M[2], M[3], M[4], b[0], b[1], b[2], b[3], b[4]

def Coefs_ydeflect(x, x_1, x_2, x_3, x_a, s, h_a, theta, E, I_zz, G, J, P, q_x, t_x, answer, d_arm):
    x_I = x_2 - (x_a / 2)
    x_II = x_2 + (x_a / 2)
    
    C_F_z1 = 0
    C_F_z2 = 0
    C_F_z3 = 0
    C_F_y1 = (-(1/6) * (1/(E * I_zz)) * (toggle(x, x_1))**3 + (1/(G * J)) * (s - (h_a/2))* d_arm * toggle(x, x_1))
    C_F_y2 = (-(1/6) * (1/(E * I_zz)) * (toggle(x, x_2))**3 + (1/(G * J)) * (s - (h_a/2))* d_arm * toggle(x, x_2))
    C_F_y3 = (-(1/6) * (1/(E * I_zz)) * (toggle(x, x_3))**3 + (1/(G * J)) * (s - (h_a/2))* d_arm * toggle(x, x_3))
    C_H = (((1/6) * (1/(E * I_zz)) * sin(theta) * (toggle(x, x_I))**3) + ((1/(G * J)) * ((h_a/2 * cos(theta)) - (s * sin(theta))) * d_arm * toggle(x, x_I)))
    C_C1 = x
    C_C2 = 1
    C_C3 = 0
    C_C4 = 0
    C_C5 = d_arm
    Coefs = [C_F_z1, C_F_z2, C_F_z3, C_F_y1, C_F_y2, C_F_y3, C_H, C_C1, C_C2, C_C3, C_C4, C_C5]
    
    b = answer - (P * (((1/6) * (1/(E * I_zz)) * sin(theta) * (toggle(x, x_II))**3) + ((1/(G * J)) * ((h_a/2 * cos(theta)) - (s * sin(theta))) * d_arm * toggle(x, x_II)))) - ((1/(E * I_zz)) * integrate(integrate(integrate(integrate(q_x,0,x),0,x),0,x),0,x)[-1]) + ((1/(G * J)) * d_arm * integrate(integrate(t_x,0,x),0,x)[-1])
    return Coefs, b

def Coefs_zdeflect(x, x_1, x_2, x_3, x_a, theta, E, I_yy, P, answer):
    x_I = x_2 - (x_a / 2)
    x_II = x_2 + (x_a / 2)
    
    C_F_z1 = ((1/6) * (1/(E * I_yy)) * (toggle(x, x_1))**3)
    C_F_z2 = ((1/6) * (1/(E * I_yy)) * (toggle(x, x_2))**3)
    C_F_z3 = ((1/6) * (1/(E * I_yy)) * (toggle(x, x_3))**3)
    C_F_y1 = 0
    C_F_y2 = 0
    C_F_y3 = 0
    C_H = (-(1/6) * (1/(E * I_yy)) * cos(theta) * (toggle(x, x_I))**3)
    C_C1 = 0
    C_C2 = 0
    C_C3 = x
    C_C4 = 1
    C_C5 = 0

    Coefs = [C_F_z1, C_F_z2, C_F_z3, C_F_y1, C_F_y2, C_F_y3, C_H, C_C1, C_C2, C_C3, C_C4, C_C5]
    
    b = answer + (1/6) * P * (1/(E * I_yy)) * cos(theta) * (toggle(x, x_II))**3
    return Coefs, b

#### TEST Inputs ####

l_a = 2.771                 #m
x_1 = 0.153                 #m
x_2 = 1.281                 #m
x_3 = 2.681                 #m
x_a = 28.0 * 10**(-2)       #m
h_a = 22.5 * 10**(-2)       #m
d_1 = 1.103 * 10**(-2)      #m
d_3 = 1.642 * 10**(-2)      #m
theta = 26 * np.pi / 180    #rad
s = 0.11427352 + h_a/2      #m

P = 91.7 * 10**3            #N

#q_x = np.linspace(1.0, 1.0, 100)     #N/m
#t_x = np.linspace(1.0, 1.0, 100)     #N
qobject = open('q.txt', 'r')
q_x = ''
for line in qobject.readlines():
    q_x += line
q_x = q_x.split(', ')
q_x[0] = q_x[0].replace('[','')
q_x[-1] = q_x[-1].replace(']','')
qobject.close()
for i in range(len(q_x)):
    q_x[i] = float(q_x[i])
q_x = np.array(q_x)
q_x = 1000 * q_x

tauobject = open('tau.txt', 'r')
t_x = ''
for line in tauobject.readlines():
    t_x += line
t_x = t_x.split(', ')
t_x[0] = t_x[0].replace('[','')
t_x[-1] = t_x[-1].replace(']','')
tauobject.close()
for i in range(len(t_x)):
    t_x[i] = float(t_x[i])
t_x = np.array(t_x)
t_x = 1000 * t_x

E = 73.1 * 10**9                   #GPa
I_zz = 1.2689896034082582 * 10**(-5)     #m^4
I_yy = 6.616599015880662 * 10**(-5)      #m^4

G = 28 * 10**9                    #GPa
J = 1.66313927 * 10**(-5)   #m^4

#### Code ####

x_I = x_2 - (x_a / 2)
x_II = x_2 + (x_a / 2)

M = np.zeros((12,12))   #matrix of equations M[row][column]
b = np.zeros(12)        #matrix of solutions of the equations
# x=[F_z1,F_z2,F_z3,F_y1,F_y2,F_y3,H,C1,C2,C3,C4,C5]

M[0], M[1], M[2], M[3], M[4], b[0], b[1], b[2], b[3], b[4] = Coefs_equi(x_1, x_2, x_3, x_a, s, h_a, l_a, theta, P, q_x, t_x)
M[5], b[5] = Coefs_ydeflect(x_1, x_1, x_2, x_3, x_a, s, h_a, theta, E, I_zz, G, J, P, q_x, t_x, (d_1 * cos(theta)), (s - h_a/2))
M[6], b[6] = Coefs_ydeflect(x_2, x_1, x_2, x_3, x_a, s, h_a, theta, E, I_zz, G, J, P, q_x, t_x, 0, (s - h_a/2))
M[7], b[7] = Coefs_ydeflect(x_3, x_1, x_2, x_3, x_a, s, h_a, theta, E, I_zz, G, J, P, q_x, t_x, (d_3 * cos(theta)), (s - h_a/2))

M[8], b[8] = Coefs_zdeflect(x_1, x_1, x_2, x_3, x_a, theta, E, I_yy, P, (-d_1 * sin(theta)))
M[9], b[9] = Coefs_zdeflect(x_I, x_1, x_2, x_3, x_a, theta, E, I_yy, P, 0)
M[10], b[10] = Coefs_zdeflect(x_2, x_1, x_2, x_3, x_a, theta, E, I_yy, P, 0)
M[11], b[11] =  Coefs_zdeflect(x_3, x_1, x_2, x_3, x_a, theta, E, I_yy, P, (-d_3 * sin(theta)))

#F= linalg.gmres(M,b,tol = 1e-5 )
F1=np.linalg.lstsq(M,b)
[F_z1,F_z2,F_z3,F_y1,F_y2,F_y3,H,C1,C2,C3,C4,C5]=F1[0]
x= np.linspace(0,l_a,1000)
My=np.zeros(len(x))
Mz=np.zeros(len(x))
T=np.zeros(len(x))
vy=np.zeros(len(x))
vz=np.zeros(len(x))
for i in range(len(x)):
   
    x1=x_1
    x2=x_2
    x3=x_3
    S=s
    My[i]= -F_z1*toggle(x[i],x1)**1-F_z2*toggle(x[i],x2)-F_z3*toggle(x[i],x3)+ H*cos(theta)*toggle(x[i],x_I)+P*cos(theta)*toggle(x[i],x_II)
   
    Mz[i]= F_y1*toggle(x[i],x1)+F_y2*toggle(x[i],x2)+F_y3*toggle(x[i],x3)
    -H*sin(theta)*toggle(x[i],x_I)-P*sin(theta)*toggle(x[i],x_II)-integrate(integrate(q_x,0,x[i]),0,x[i])[-1]
    
    T[i]=(F_y1*(S-(h_a/2))*toggle0(x[i],x1)+F_y2*(S-(h_a/2))*toggle0(x[i],x2)+F_y3*(S-(h_a/2))*toggle0(x[i],x3)
   -P*sin(theta)*S*toggle0(x[i],x_II)-H*sin(theta)*S*toggle0(x[i],x_I)+P*cos(theta)*(h_a/2)*toggle0(x[i],x_II)+H*cos(theta)*(h_a/2)*toggle0(x[i],x_I)-integrate(t_x,0,x[i]))[-1]*-1
    vy[i]=-(1/(E*I_zz))*((1/6)*F_y1*toggle(x[i],x1)**3 +(1/6)*F_y2*toggle(x[i],x2)**3+(1/6)*F_y3*toggle(x[i],x3)**3-(1/6)*H*sin(theta)*toggle(x[i],x_I)**3-(1/6)*P*sin(theta)*toggle(x[i],x_II)**3-integrate(integrate(integrate(integrate(q_x,0,x[i]),0,x[i]),0,x[i]),0,x[i])[-1])+C1*x[i]+C2
    vz[i]=-(1/(E*I_yy))*(-(1/6)*F_z1*toggle(x[i],x1)**3-(1/6)*F_z2*toggle(x[i],x2)**3-(1/6)*F_z3*toggle(x[i],x3)**3+(1/6)*H*cos(theta)*toggle(x[i],x_I)**3+(1/6)*P*cos(theta)*toggle(x[i],x_II)**3)+C3*x[i]+C4

d=integrate(integrate(My,0,l_a)+C1,0,l_a) 
d=np.array(d)
d=d*(1/(E*I_zz)) 
#plt.plot(,'green')
#plt.plot(x[0:len(d)],d+C2,'green')
plt.plot(x[0:len(vz)],vz,'red')
#plt.plot(x,np.zeros(len(x)))

plt.show()