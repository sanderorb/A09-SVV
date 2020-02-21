import numpy as np
import matplotlib.pyplot as plt
import math as m
from mpl_toolkits import mplot3d
import scipy.integrate as integrate
#==========Inputs==========#
#Variables
Ca=0.547 #[m]
la=2.771 #[m]
Nz=81
Nx=41

zsc=-2 #Shear center z coordinate

#Aero Load
data = np.genfromtxt("qload.txt", delimiter=",")

#New grid for new aero load, moments, shear at new positions
Nxnew=400
Nznew=60 #Nz only needed for a aero load with new grid
#Nznew > Nxnew usually

#==========Extracting coordinates==========#
thetaz=[]
thetax=[]
zcoord=[]
xcoord=[]

for i in range(1,Nz+2):
    t=(i-1)/Nz*m.pi
    thetaz.append(t)
for i in range(0,len(thetaz)-1):
    z=-0.5*(Ca/2*(1-m.cos(thetaz[i]))+Ca/2*(1-m.cos(thetaz[i+1])))
    zcoord.append(z)

for i in range(1,Nx+2):
    t=(i-1)/Nx*m.pi
    thetax.append(t)
for i in range(0,len(thetax)-1):
    x=0.5*(la/2*(1-m.cos(thetax[i]))+la/2*(1-m.cos(thetax[i+1])))
    xcoord.append(x)

#==========Patch Interpolation==========#
#Unit square coordinates matrix
A=np.array([[1,-1, 1,-1],
            [1, 1, 1, 1],
            [1,-1,-1, 1],
            [1, 1,-1,-1]])

coefs = np.zeros((80, 40, 4)) #Empty matrix of coefficients (80x40 segments not points)
for i in range(0,80):
    for j in range(0,40):
        f=np.array([[data[i,j]],[data[i,j+1]],[data[i+1,j]],[data[i+1,j+1]]])
        M=np.matmul(np.linalg.inv(A),f) #Returns coefficients a0, a1, a2 and a3 from a0+a1*x+a2*y+a3*xy
        M=np.transpose(M)
        mx=2/(xcoord[j+1]-xcoord[j]) #Scaling of coordinates
        cx=1-mx*xcoord[j+1]
        mz=2/(zcoord[i]-zcoord[i+1])
        cz=1-mz*zcoord[i]
        M[0][0]=M[0][0]+M[0][1]*cx+M[0][2]*cz+M[0][3]*cx*cz #Applied scaling to polynomial coefficients
        M[0][1]=M[0][1]*mx+M[0][3]*mx*cz
        M[0][2]=M[0][2]*mz+M[0][3]*mz*cx
        M[0][3]=M[0][3]*mx*mz
        coefs[i,j]=M

#==========Aerodynamic Load with new grid==========#
xcoordnew=np.linspace(xcoord[0],xcoord[40],Nxnew)
zcoordnew=np.linspace(zcoord[0],zcoord[80],Nznew)
qnew = np.zeros((Nznew, Nxnew))
for j in range(0,Nznew):
    for i in range(0,Nxnew):
        indx=np.argmin(abs(xcoord-xcoordnew[i]))-1
        indz=np.argmin(abs(zcoord-zcoordnew[j]))-1
        if indx<0:
            indx=0
        if indz<0:
            indz=0
        if indx==41:
            indx=40
        if indz==81:
            indz=80
        c0, c1, c2, c3=coefs[indz,indx]
        qnew[j,i]=c0+c1*xcoordnew[i]+c2*zcoordnew[j]+c3*xcoordnew[i]*zcoordnew[j]

#==========Intergration in z for w(x)==========#
#Obtaining new coefficients
wcoefs = np.zeros((80, 40, 2))
for i in range(0,40):
    for j in range(0,80):
        a=coefs[j,i]
        b1=a[0]*zcoord[j+1]+0.5*a[2]*zcoord[j+1]**2-(a[0]*zcoord[j]+0.5*a[2]*zcoord[j]**2)
        b2=a[1]*zcoord[j+1]+0.5*a[3]*zcoord[j+1]**2-(a[1]*zcoord[j]+0.5*a[3]*zcoord[j]**2)
        b=[b1,b2] #Results in polynomial F(x)=b1+b2*x
        wcoefs[j,i]=b

print(wcoefs)
#Summing over one chord
wresults= []
for i in range(40):
    columns= wcoefs[:,i]
    wresults.append([sum(columns[:,0]), sum(columns[:,1])])
    if i+1==40:
        columns= wcoefs[:,i]
        wresults.append([sum(columns[:,0]), sum(columns[:,1])])
print(wresults)
#Calculating values
wresultsnew=[]
for i in range(0,Nxnew):
    indx=np.argmin(abs(xcoord-xcoordnew[i]))-1
    if indx<0:
        indx=0
    if indx==41:
        indx=40
    else:
        wresultsnew.append(wresults[indx][0]+wresults[indx][1]*xcoordnew[i])

#==========Intergration in z for tau(x)==========#
taucoefs = np.zeros((80, 40, 2))
for i in range(0,40):
    for j in range(0,80):
        a=coefs[j,i]
        t0=0.5*a[0]*zcoord[j]**2+a[2]/3*zcoord[j]**3-zsc*(a[0]*zcoord[j]+0.5*a[2]*zcoord[j]**2)
        t1=0.5*a[0]*zcoord[j+1]**2+a[2]/3*zcoord[j+1]**3-zsc*(a[0]*zcoord[j+1]+0.5*a[2]*zcoord[j+1]**2)
        t2=0.5*a[1]*zcoord[j]**2+a[3]/3*zcoord[j]**3-zsc*(a[1]*zcoord[j]+0.5*a[3]*zcoord[j]**2)
        t3=0.5*a[1]*zcoord[j+1]**2+a[3]/3*zcoord[j+1]**3-zsc*(a[1]*zcoord[j+1]+0.5*a[3]*zcoord[j+1]**2)
        b=[t1-t0,t3-t2] #Results in polynomial F(x)=b1+b2*x
        taucoefs[j,i]=b



tauresults= []
for i in range(40):
    columns= taucoefs[:,i]
    tauresults.append([sum(columns[:,0]), sum(columns[:,1])])
    if i+1==40:
        columns= taucoefs[:,i]
        tauresults.append([sum(columns[:,0]), sum(columns[:,1])])

tauresultsnew=[]
for i in range(0,Nxnew):
    indx=np.argmin(abs(xcoord-xcoordnew[i]))-1
    if indx<0:
        indx=0
    if indx==41:
        indx=41
    else:
        tauresultsnew.append(tauresults[indx][0]+tauresults[indx][1]*xcoordnew[i])

# print(wresultsnew)
# print(tauresultsnew)

# #==========Plotting==========#
#Grids for plotting
xcoordgrid, zcoordgrid = np.meshgrid(xcoord, zcoord)
xcoordgridnew, zcoordgridnew = np.meshgrid(zcoordnew, xcoordnew)

#3D plot
#ax = plt.axes(projection='3d')

#Grid plot
#plt.scatter(xcoordgrid, zcoordgrid)
#plt.scatter(xcoordgridnew, zcoordgridnew)

#Aero load plots
#ax.plot_surface(xcoordgrid, zcoordgrid, data, rstride=1, cstride=1,cmap='viridis', edgecolor='none') #Lift distribution
#ax.plot_surface(xcoordgridnew, zcoordgridnew, np.transpose(qnew), rstride=1, cstride=1,cmap='viridis', edgecolor='none') #Lift distribution

#Force and Moment plots
#plt.scatter(xcoordnew, wresultsnew)
plt.scatter(xcoordnew, tauresultsnew)

plt.show()
