import matplotlib.pyplot as plt
import numpy as np

l_a = 2.771 #m

Myobject = open('My.txt', 'r')
My = ''
for line in Myobject.readlines():
    My += line
My = My.replace('[','')
My = My.replace(']','')
My = My.split(' ')
My = list(filter(('').__ne__, My))
for i in range(len(My)):
    My[i] = float(My[i])
Myobject.close()

Mzobject = open('Mz.txt', 'r')
Mz = ''
for line in Mzobject.readlines():
    Mz += line
Mz = Mz.replace('[','')
Mz = Mz.replace(']','')
Mz = Mz.split(' ')
Mz = list(filter(('').__ne__, Mz))
for i in range(len(Mz)):
    Mz[i] = float(Mz[i])
Mzobject.close()

Syobject = open('Sy.txt', 'r')
Sy = ''
for line in Syobject.readlines():
    Sy += line
Sy = Sy.replace('[','')
Sy = Sy.replace(']','')
Sy = Sy.split(' ')
Sy = list(filter(('').__ne__, Sy))
for i in range(len(Sy)):
    Sy[i] = float(Sy[i])
Syobject.close()

Szobject = open('Sz.txt', 'r')
Sz = ''
for line in Szobject.readlines():
    Sz += line
Sz = Sz.replace('[','')
Sz = Sz.replace(']','')
Sz = Sz.split(' ')
Sz = list(filter(('').__ne__, Sz))
for i in range(len(Sz)):
    Sz[i] = float(Sz[i])
Szobject.close()

Tobject = open('T.txt', 'r')
T = ''
for line in Tobject.readlines():
    T += line
T = T.replace('[','')
T = T.replace(']','')
T = T.split(' ')
T = list(filter(('').__ne__, T))
for i in range(len(T)):
    T[i] = float(T[i])
Tobject.close()

vyobject = open('vy.txt', 'r')
vy = ''
for line in vyobject.readlines():
    vy += line
vy = vy.replace('[','')
vy = vy.replace(']','')
vy = vy.split(' ')
vy = list(filter(('').__ne__, vy))
for i in range(len(vy)):
    vy[i] = float(vy[i])
vyobject.close()

vzobject = open('vz.txt', 'r')
vz = ''
for line in vzobject.readlines():
    vz += line
vz = vz.replace('[','')
vz = vz.replace(']','')
vz = vz.split(' ')
vz = list(filter(('').__ne__, vz))
for i in range(len(vz)):
    vz[i] = float(vz[i])
vzobject.close()

phiobject = open('phi.txt', 'r')
phi = ''
for line in phiobject.readlines():
    phi += line
phi = phi.replace('[','')
phi = phi.replace(']','')
phi = phi.split(' ')
phi = list(filter(('').__ne__, phi))
for i in range(len(phi)):
    phi[i] = float(phi[i])
phi = list(np.array(phi) * 180/np.pi * -1)
phiobject.close()

Mzmodelobject = open('Mzmodel.txt', 'r')
Mzmodel = ''
for line in Mzmodelobject.readlines():
    Mzmodel += line
Mzmodel = Mzmodel.replace('[','')
Mzmodel = Mzmodel.replace(']','')
Mzmodel = Mzmodel.split(', ')
Mzmodel = list(filter(('').__ne__, Mzmodel))
for i in range(len(Mzmodel)):
    Mzmodel[i] = float(Mzmodel[i])
Mzmodel = list(reversed(Mzmodel))
Mzmodelobject.close()

Mymodelobject = open('Mymodel.txt', 'r')
Mymodel = ''
for line in Mymodelobject.readlines():
    Mymodel += line
Mymodel = Mymodel.replace('[','')
Mymodel = Mymodel.replace(']','')
Mymodel = Mymodel.split(', ')
Mymodel = list(filter(('').__ne__, Mymodel))
for i in range(len(Mymodel)):
    Mymodel[i] = float(Mymodel[i])
Mymodelobject.close()

Szmodelobject = open('Szmodel.txt', 'r')
Szmodel = ''
for line in Szmodelobject.readlines():
    Szmodel += line
Szmodel = Szmodel.replace('[','')
Szmodel = Szmodel.replace(']','')
Szmodel = Szmodel.split(', ')
Szmodel = list(filter(('').__ne__, Szmodel))
for i in range(len(Szmodel)):
    Szmodel[i] = float(Szmodel[i])
Szmodelobject.close()

Symodelobject = open('Symodel.txt', 'r')
Symodel = ''
for line in Symodelobject.readlines():
    Symodel += line
Symodel = Symodel.replace('[','')
Symodel = Symodel.replace(']','')
Symodel = Symodel.split(', ')
Symodel = list(filter(('').__ne__, Symodel))
for i in range(len(Symodel)):
    Symodel[i] = float(Symodel[i])
Symodelobject.close()

Tmodelobject = open('Tmodel.txt', 'r')
Tmodel = ''
for line in Tmodelobject.readlines():
    Tmodel += line
Tmodel = Tmodel.replace('[','')
Tmodel = Tmodel.replace(']','')
Tmodel = Tmodel.split(', ')
Tmodel = list(filter(('').__ne__, Tmodel))
for i in range(len(Tmodel)):
    Tmodel[i] = float(Tmodel[i])
Tmodelobject.close()

vymodelobject = open('vymodel.txt', 'r')
vymodel = ''
for line in vymodelobject.readlines():
    vymodel += line
vymodel = vymodel.replace('[','')
vymodel = vymodel.replace(']','')
vymodel = vymodel.split(', ')
vymodel = list(filter(('').__ne__, vymodel))
for i in range(len(vymodel)):
    vymodel[i] = float(vymodel[i])
vymodelobject.close()

vzmodelobject = open('vzmodel.txt', 'r')
vzmodel = ''
for line in vzmodelobject.readlines():
    vzmodel += line
vzmodel = vzmodel.replace('[','')
vzmodel = vzmodel.replace(']','')
vzmodel = vzmodel.split(', ')
vzmodel = list(filter(('').__ne__, vzmodel))
for i in range(len(vzmodel)):
    vzmodel[i] = float(vzmodel[i])
vzmodelobject.close()

phimodelobject = open('phimodel.txt', 'r')
phimodel = ''
for line in phimodelobject.readlines():
    phimodel += line
phimodel = phimodel.replace('[','')
phimodel = phimodel.replace(']','')
phimodel = phimodel.split(', ')
phimodel = list(filter(('').__ne__, phimodel))
for i in range(len(phimodel)):
    phimodel[i] = float(phimodel[i])
phimodel = list(np.array(phimodel) * 180/np.pi)
phimodelobject.close()

####

xtab = np.linspace(0, l_a, 1000)

plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)
plt.rc('axes', labelsize=15)
plt.rc('figure', titlesize=40)
plt.plot(xtab, My, 'red')
plt.plot(xtab, Mymodel, 'blue')
plt.title('My')
plt.xlabel('m')
plt.ylabel('Nm')
plt.show()

plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)
plt.rc('axes', labelsize=15)
plt.rc('figure', titlesize=40)
plt.plot(xtab, Mz, 'red')
plt.plot(xtab, Mzmodel, 'blue')
plt.title('Mz')
plt.xlabel('m')
plt.ylabel('Nm')
plt.show()

plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)
plt.rc('axes', labelsize=15)
plt.rc('figure', titlesize=40)
plt.plot(xtab, Sy, 'red')
plt.plot(xtab, Symodel, 'blue')
plt.title('Sy')
plt.xlabel('m')
plt.ylabel('N')
plt.show()

plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)
plt.rc('axes', labelsize=15)
plt.rc('figure', titlesize=40)
plt.plot(xtab, Sz, 'red')
plt.plot(xtab, Szmodel, 'blue')
plt.title('Sz')
plt.xlabel('m')
plt.ylabel('N')
plt.show()

plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)
plt.rc('axes', labelsize=15)
plt.rc('figure', titlesize=40)
plt.plot(xtab, T, 'red')
plt.plot(xtab, Tmodel, 'blue')
plt.title('T')
plt.xlabel('m')
plt.ylabel('Nm')
plt.show()

plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)
plt.rc('axes', labelsize=15)
plt.rc('figure', titlesize=40)
plt.plot(xtab, vy, 'red')
plt.plot(xtab, vymodel, 'blue')
plt.title('vy')
plt.xlabel('m')
plt.ylabel('m')
plt.show()

plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)
plt.rc('axes', labelsize=15)
plt.rc('figure', titlesize=40)
plt.plot(xtab, vz, 'red')
plt.plot(xtab, vzmodel, 'blue')
plt.title('vz')
plt.xlabel('m')
plt.ylabel('m')
plt.show()

plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)
plt.rc('axes', labelsize=15)
plt.rc('figure', titlesize=40)
plt.plot(xtab, phi, 'red')
plt.plot(xtab, phimodel, 'blue')
plt.title('phi')
plt.xlabel('m')
plt.ylabel('deg')
plt.show()
