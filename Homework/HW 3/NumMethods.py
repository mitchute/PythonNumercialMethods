"""
Numerical Methods
"""

import numpy as np

###########################################

def lagrange(xd, fd, x):
    f = 0.0
    for j in range(len(xd)):
        b = 1.0
        for i in range(len(xd)):
            if i != j:
                b = b * ((x - xd[i])/(xd[j] - xd[i]))
        f = f + fd[j] * b
    return f

###########################################

def lagrange_derivative(xd, fd, x):
    a = np.zeros(len(xd))
    b = np.zeros(len(xd))
    d = np.zeros(len(xd))

    #Precompute the denominator
    for j in range(len(xd)):
        a[j] = 1.0
        for i in range(len(xd)):
            if i != j:
                a[j] = a[j] * (xd[j] - xd[i])

    #d(j) is the base function for the derivative
    #derivative of the lagrange polynomial
    p = 0.0
    for j in range(len(xd)):
        d[j] = 0.0
        for i in range(len(xd)):
            b[i] = 1.0
            if i != j:
                for k in range(len(xd)):
                    if k != i and k != j:
                        b[i] = b[i] * (x - xd[k])
                    # end if
                # end for
                d[j] = d[j] + b[i]/a[j]
            # end if
        # end for
        p = p + d[j] * fd[j]
    # end for

    return p

###########################################

def c2(x, u):

    n = len(u)

    h = x[1]-x[0]

    up = np.zeros(n)

    for i in range(1, (n - 1)):
        up[i] = (u[i+1] - u[i-1])/(2.0*h)

    # sided difference for i == 1
    i = 0
    up[i] = (-3.0*u[i] + 4.0*u[i+1] - u[i+2])/(2.0*h)

    # sided difference for i == n
    i = n - 1
    up[i] = (-3.0*u[i] + 4.0*u[i-1] - u[i-2])/(-2.0*h)

    return up

###########################################

def pade4_twice(x, u):

    # calls pade4 twice to approximate the second derivative
    d1 = pade4(x, u)

    return pade4(x, d1)

###########################################

def pade4_2(x, u):

    # 4th order compact Pade scheme for the second order derivative(up)
    # 3-4-3
    # 3rd order b.c.
    # 4th order interior

    n = len(u)

    h = x[1]-x[0]

    a = np.zeros(n)
    b = np.zeros(n)
    c = np.zeros(n)
    r = np.zeros(n)

    i = 0
    b[i] = 1.0
    c[i] = 11.0
    r[i] = (13.0*u[i] - 27.0*u[i+1] + 15*u[i+2] - u[i+3])/(h**2)

    for i in range(1, (n - 1)):
        a[i] = 1.0/12.0
        b[i] = 10.0/12.0
        c[i] = 1.0/12.0
        r[i] = (1.0/h**2)*(u[i-1] - 2*u[i] + u[i+1])

    i = n - 1
    a[i] = 11.0
    b[i] = 1.0
    r[i] = (13.0*u[i] - 27.0*u[i-1] + 15*u[i-2] - u[i-3])/(h**2)

    return tdma(a, b, c, r)

###########################################

def pade4(x, u):

    # 4th order compact Pade scheme for the first order derivative(up)
    # 3-4-3
    # 3rd order b.c.
    # 4th order interior

    n = len(u)

    h = x[1]-x[0]

    a = np.zeros(n)
    b = np.zeros(n)
    c = np.zeros(n)
    r = np.zeros(n)

    i = 0
    b[i] = 1.0
    c[i] = 2.0
    r[i] = (-5.0*u[i] + 4.0*u[i+1] + u[i+2])/(2.0*h)

    for i in range(1, (n - 1)):
        a[i] = 1.0
        b[i] = 4.0
        c[i] = 1.0
        r[i] = (3.0/h)*(u[i+1] - u[i-1])

    i = n - 1
    a[i] = 2.0
    b[i] = 1.0
    r[i] = (-5.0*u[i] + 4.0*u[i-1] + u[i-2])/(-2.0*h)

    return tdma(a, b, c, r)

###########################################

def tdma(a, b, c, r):

    # Tridiagonal matrix algorithm (TDMA)
    # Thomas algorithm
    # Solution tridiagonal systems
    # a: lower diagonal
    # b: Main diagonal
    # c: upper diagonal
    # r: source vector
    # x: solution vector
    #   for indices s(start) to e(end)
    #   i: s,s+1,s+2,...,i,...,e
    #
    # Note: a(s) and c(e) are dummy coefficients, not used.

    n = len(r)
    s = 0
    e = n - 1

    x = np.zeros(n)

    # Forward elmination phase
    for i in range(1, n):
        b[i] = b[i] - a[i]/b[i-1]*c[i-1]
        r[i] = r[i] - a[i]/b[i-1]*r[i-1]

    # Backward substitution phase
    x[e] = r[e]/b[e]
    for i in range((e-1), (s-1), -1):
        x[i] = (r[i] - c[i]*x[i+1])/b[i]

    return x

###########################################

def trap1D(x, f):

    # Trapezoidal ruld for numerical integration of f(x)
    # for equally spaced distributed mesh with interval dx
    # Second-order accurate

    dx = x[1]-x[0]
    th = dx/2.0
    n = len(x)

    sum = 0.0
    for i in range(n-1):
        ds = th * (f[i]+f[i+1])
        sum = sum + ds

    return sum

###########################################

def trapEnd1D(x, f):

    # Trapezoidal ruld for numerical integration of f(x)
    # plus end correction
    # for equally spaced distributed mesh with interval dx
    # Second-order accurate

    dx = x[1] - x[0]
    sum = trap1D(x, f)

    # end correction
    f_der = pade4(x, f)
    sum = sum - dx**2/12.0*(f[-1]-f[0])

    return sum

###########################################

def simp1D(x, f):

    # Simpson's 1/3 rul for numerical intgration of f(x)
    # for equally distributed mesh with interval dx
    # nx should even number
    # fourth-order accurate

    dx = x[1]-x[0]
    n = len(x)
    nh = int(n/2)
    th = 1.0/3.0*dx

    sum = 0.0
    for i in range(nh):
        ds = th * (f[2*i] + 4.0*f[2*i+1] + f[2*i+2])
        sum = sum + ds

    return sum

###########################################

def gaussQuadrature(a, b):

    def f(x):
        return np.math.exp(x)

    k = []
    w = []

    k.append(0.0)
    k.append(1.0/3.0 * np.math.sqrt(5.0 - 2.0*np.math.sqrt(10.0/7.0)))
    k.append(-1.0/3.0 * np.math.sqrt(5.0 - 2.0*np.math.sqrt(10.0/7.0)))
    k.append(1.0/3.0 * np.math.sqrt(5.0 + 2.0*np.math.sqrt(10.0/7.0)))
    k.append(-1.0/3.0 * np.math.sqrt(5.0 + 2.0*np.math.sqrt(10.0/7.0)))

    w.append(128.0/225.0)
    w.append((322.0 + 13.0*np.math.sqrt(70.0))/900.0)
    w.append((322.0 + 13.0*np.math.sqrt(70.0))/900.0)
    w.append((322.0 - 13.0*np.math.sqrt(70.0))/900.0)
    w.append((322.0 - 13.0*np.math.sqrt(70.0))/900.0)

    sum = 0.0

    for i in range(4):
        sum = sum + (b-a)/2.0*f((b+a)/2.0 + (b-a)/2.0*k[i])*w[i]

    return sum

###########################################

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
""" Right hand side (RHS) of equations """

def RHS(x,y):
    # Right hand side
    n = len(y)
    r = np.zeros(n)

    a = -(x + 3)/(x + 1)
    b = (x + 3)/(x + 1)**2
    f = 2*(x + 1) +3*b

    r[0] = y[1]
    r[1] = f - a*y[1] - b*y[0]

    return r

###########################################

def print_hello():
	return 'Hello Numerical Methods User!'

###########################################

if __name__ == "__main__":
	# Unit test type stuff here
    print( print_hello() + "\n")

    f_x = [0,1,2,2,2,1,0]
    x = [0,1,2,3,4,5,6]

    print("2nd Order Central Diff Derrivative for f(x) =: "
         + str(f_x) + " with a uniform step size is: "
         + str(c2(x, f_x)) + "\n")

    print("4nd Order Pade Scheme 1st-order Derivative for f(x) =: "
         + str(f_x) + " with a uniform step size is: "
         + str(pade4(x, f_x)) + "\n")

    print("4nd Order Pade Scheme 2nd-order Derivative for f(x) =: "
         + str(f_x) + " with a uniform step size is: "
         + str(pade4_2(x, f_x)) + "\n" )

    print("1D Trapezoidal integration for f(x) =: "
         + str(f_x) + " and x =: " + str(x) + " is: "
         + str(trap1D(x, f_x)) + "\n" )

    print("1D Trapezoidal integration with end correction for f(x) =: "
         + str(f_x) + " and x =: " + str(x) + " is: "
         + str(trapEnd1D(x, f_x)) + "\n" )

    print("1/3 Simpson's rule integration for f(x) =: "
         + str(f_x) + " and x =: " + str(x) + " is: "
         + str(simp1D(x, f_x)) + "\n" )
