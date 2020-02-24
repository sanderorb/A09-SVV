import numpy as np
import math as m

Ca = 0.547  # m
ha = 0.225  # m
r = ha/2.
h = r
tsk = 1.1/1000.  # m
tsp = 2.9/1000.  # m

diag = m.sqrt((Ca-r)**2+h**2)

qb1 = -0.0003439746441695739*10**3
qb2 = -0.0002772647948510785*10**3
qb3 = -0.00084389730478562*10**3


Am1 = (1/2)*m.pi*r**2
Am2 = h*(Ca-r)

c1 = 1/(2*Am1)
c2 = 1/(2*Am2)

a1 = c1*((m.pi*r/tsk)+(ha/tsp))
a2 = c1*(-ha/tsp)
b1 = c2*(-ha/tsp)
b2 = c2*((diag*2/tsp)+(ha/tsp))

k1 = -c1*((m.pi*r*qb1/tsk)-(ha*qb2/tsp))
k2 = -c2*((2*diag*qb3/tsk)+(ha*qb2/tsp))

A = np.array([[a1,a2],[b1,b2]])
b = np.array([[k1],[k2]])
z = np.linalg.solve(A,b)

print("[q01,q02]= ",z[0],z[1])

q1 = z[0]+2*qb1         #Shear flow LE skin
q2 = z[0]-z[1]-2*qb2    #Shear flow Spar
q3 = z[1]+qb3           #Shear flow TE skin

Sy = 1                  #Unit Load
arm = (h*diag)/(Ca-r)     #Arm for, Shear flow TE skin

eta = (((q1*m.pi*r*0.5))*r+((q3*diag)*2*arm))/Sy    #Moment Equivalent

print("eta=", eta)


#------------------------------------------
#Torsional stiffness constant
d1 = 2*Am1
d2 = 2*Am2
d3 = 0
e1 = c1*((m.pi*r/tsk)+ha/tsp)   #*q0.1
e2 = c1*(-ha/tsp)               #*q0.2
e3 = -1                         #*Gdtheta/dz
f1 = c2*(-ha/tsp)               #*q0.1
f2 = c2*((2*diag/tsk)+ha/tsp)   #*q0.2
f3 = -1                         #*Gdtheta/dz

T = 1

B = np.array([[d1,d2,d3],[e1,e2,e3],[f1,f2,f3]])
c = np.array([[T],[0],[0]])

q = np.linalg.solve(B,c)        #[q0.1,q0.2,Gdtheta/dz]

J = T*1000/q[2]
print()
print("shear flow=",q[0],q[1])
print("Torsional Stiffness=", J)
