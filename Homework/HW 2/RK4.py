import numpy as np

def RK4(x,h,y):
    # Runga-Kutta fourth-order scheme

    n = len(y)
    k1 = np.zeros(len(y))
    k2 = np.zeros(len(y))
    k3 = np.zeros(len(y))
    k4 = np.zeros(len(y))
    r = np.zeros(len(y))

    r = RHS(x,y)

    for k in range(n):
        k1[k] = h*r[k]

    r = RHS(x+h/2.0,y+k1/2.0)

    for k in range(n):
        k2[k] = h*r[k]

    r = RHS(x+h/2.0,y+k2/2.0)

    for k in range(n):
        k3[k] = h*r[k]

    r = RHS(x+h, y+k3)

    for k in range(n):
        k4[k] = h*r[k]

    for k in range(n):
        y[k] = y[k] + (k1[k] + 2.0*(k2[k] + k3[k]) + k4[k])/6.0

    return y

###########################################

def RHS(x,y):
    # Right hand side of Blasius Equation
    n = len(y)
    r = np.zeros(n)

    r[0] = -y[0]*y[2]
    r[1] = y[0]
    r[2] = y[1]

    return r