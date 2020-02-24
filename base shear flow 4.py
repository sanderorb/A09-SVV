    # Base shear flow calculation
import math as m
import numpy as np


    # Moments of inertia
Izz = 1 
Iyy = 1
strspac = 73.5932358            # stringer spacing
t = 0.1                         # skin thickness
ts = 0.5                        # spar thickness
h = 5                           # half length of spar
l = 2.771 * 10**3
    # Boom properties
yb = np.array([[0.0, 0.0],
          [23.22471754086723, 68.45563484361286],
          [83.30975896608109, 108.64704242813845],
          [155.15843148890556, 101.45495157076668],
          [226.4023530363773, 83.0085967397182],
          [297.646274583849, 64.5622419086697],
          [368.8901961313207, 46.11588707762121],
          [440.1341176787924, 27.669532246572714],
          [511.37803922626415, 9.223177415524233],
          [23.22471754086723, -68.45563484361286],
          [83.30975896608109, -108.64704242813845],
          [155.15843148890556, -101.45495157076668],
          [226.4023530363773, -83.0085967397182],
          [297.646274583849, -64.5622419086697],
          [368.8901961313207, -46.11588707762121],
          [440.1341176787924, -27.669532246572714],
          [511.37803922626415, -9.223177415524233]])
b = 2 # boom area 
qb = [] # empty list to append base shear flow
d = []
# Skin contribution functions (vertical)
def s1(x):
    f = t*h*m.sin(x/h) 
    return f
def s2(x):
    f = ts*x
    return f
def s3(x):
    f = t*(h-x*(h/l))
    return f
def s4(x):
    f = t*(-x*(h/l))
    return f
def s5(x):
    f = ts*(x-h)
    return f
def s6(x):
    f = t*h*m.sin(x/h)
    return f


# Skin contribution functions (horizontal)
r = h
centroid = ([[0,213.50105304149162/1000.]])
l = Ca-r
xc = centroid[1]
diag = m.sqrt((Ca-r)**2+h**2)

def x1(y):
    f1 = t*(-(xc-r)+(r)*m.cos(y/r))
    return f1
def x2(y):
    f1 = ts*(-(xc-r))
    return f1
def x3(y):
    f1 = t*(-(xc-r)+y*(Ca-r)/diag)
    return f1
def x4(y):
    f1 = t*(-(xc-r)+(Ca-r)-y*(Ca-r)/Diag)
    return f1
def x5(y):
    f1 = ts*(-(xc-r))
    return f1
def x6(y):
    t*(-(xc-r)+(-r)*m.cos(y/r))
    return f1

# Integration function

n = 10 

def integrate(f,n,a,b):
    # f = function
    # n = number of trapezoids
    # a = startpoint
    # b = endpoint
    h = (b-a)/ float(n)
    inte = 0.5 * h * (f(0) + f(1))
    for i in range(1, int(n)):
        inte = inte + h * f(a + i*h)
    return inte

# Total skin contribution per section

S1 = integrate(s1, n, 0, h*m.pi*0.5)
S2 = integrate(s2, n, 0, h)
S3 = integrate(s3, n, 0, l)
S4 = integrate(s4, n, 0, l)
S5 = integrate(s5, n, 0, h)
S6 = integrate(s6, n, 0, h*m.pi*0.5)

SH1 = integrate(x1, n, 0, m.pi*r*0.5)
SH2 = integrate(x2, n, 0, h)
SH3 = integrate(x3, n, 0, diag)
SH4 = integrate(x1, n, 0, diag)
SH5 = integrate(x2, n, 0, h)
SH6 = integrate(x3, n, 0, m.pi*r*0.5)


    # Base shear flow calculations
for s in range(1252):
    # Boom contribution
    o = int(s/strspac)  #number of booms passed
    B = sum(b*yb[0,2][0])
    
# Skin contribution
    
    # Outer circumference
        # section 1
    if s <= 112.5*m.pi: 
        q = integrate(s1, n, 0, s/(2*h))
        # section 3
    elif s <= 625.5:   
        q = S1 + S2 + integrate(s3, n, 0, s - 112.5*m.pi)
        # section 4
    elif s <= 897.57:
        q = S1 + S2 + S3 + integrate(s4, n, 0, s - 897.57)
        # section 6
    elif s <= 1251:
        q = S1 + S3 + S4 + integrate(s6, n, 0.5*m.pi/h, (s-897.57)/h)
    else:
        print('Boundary error')
    
       
    # Total base shear flow circumference
    qt = -1/Iyy * (q + B)
    qb.append(qt)
    d.append(s)

# Spar section
    qs = []
    p = []
    for ss in range(h):
        # section 2
        qs1 = integrate(s2, n, 0, ss) 
        qs.append(qs1)
        p.append(ss)
        
print('qb circumference:', qb)
