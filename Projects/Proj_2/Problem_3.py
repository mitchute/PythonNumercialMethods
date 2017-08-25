import numpy as np
import matplotlib.pyplot as plt

################################################################################

Pi = np.math.pi
alpha = 1 / Pi**2.0

################################################################################

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

    print(a)
    print(b)
    print(c)
    print(x)

    return x

################################################################################

def exact(t,x):
    return -np.math.sin(Pi*x) * np.math.exp(-t)

################################################################################

def calcError(t,x,u):

    r = np.zeros(len(u))

    for i in range(len(u)):

        ue = exact(t,x[i])
        r[i] = abs(ue - u[i])
    return r

################################################################################

def FTCS(dt, dx, u):

    r = np.zeros(len(u))

    for i in range(len(u)):
        if i == 0:
            r[i] = 0
        elif i == len(u) - 1:
            r[i] = 0
        else:
            r[i] = u[i] + ((dt * alpha)/dx**2.0) * (u[i+1] - 2*u[i] + u[i-1])
    return r

################################################################################

def BTCS(dt, dx, u):

    # Initialize tri-diagonal system
    a = np.zeros(len(u))
    b = np.zeros(len(u))
    c = np.zeros(len(u))
    r = np.zeros(len(u))

    # Calculate Coefficients
    LowerDiag = -(alpha*dt)/(dx**2.0)
    MainDiag = (1 + 2*alpha*dt/dx**2.0)
    UpperDiag = -(alpha*dt)/(dx**2.0)

    # Set Coefficient Matrix
    for i in range(len(u) - 1):
        a[i] = LowerDiag
        b[i] = MainDiag
        c[i] = UpperDiag
        r[i] = u[i]       

    # Set Boundary Conditions
    r[0] = 0 # Left side
    r[-1] = 0 # Right side

    return tdma(a,b,c,r)

################################################################################

def CN(dt, dx, u):
    # Initialize tri-diagonal system
    a = np.zeros(len(u))
    b = np.zeros(len(u))
    c = np.zeros(len(u))
    r = np.zeros(len(u))

    # Calculate Coefficients
    LowerDiag = -(alpha*dt)/(2.0*dx**2.0)
    MainDiag = (1 + alpha*dt/dx**2.0)
    UpperDiag = -(alpha*dt)/(2.0*dx**2.0)
    Source = (alpha*dt)/(2.0*dx**2.0)

    # Set Coefficient Matrix
    a.fill(LowerDiag)
    b.fill(MainDiag)
    c.fill(UpperDiag)

    # Set Source Matrix
    for i in range(1, len(u)-1):
        r[i] = Source * (u[i+1] - 2*u[i] + u[i-1]) + u[i]

    # Set Boundary Conditions
    r[0] = 0 # Left side
    r[-1] = 0 # Right side

    return tdma(a,b,c,r)

################################################################################

def main(case):
    # Select case
    if case == 1:
        dt = 0.01
        dx = 0.05
    elif case == 2:
        dt = 0.1
        dx = 0.1
    elif case == 3:
        dt = 0.1
        dx = 0.05
    elif case == 4:
        dt = 0.001
        dx = 0.02
    elif case == 5:
        dt = 0.001
        dx = 0.025

    a = -1
    b = 1

    t = 0.0
    tol = 0.00005

    n = (b - a)/dx + 1

    # Initialize
    x = np.zeros(n)
    u = np.zeros(n)
    u_init = np.zeros(n)
    u_exact = np.zeros(n)

    # Set initial x and u
    for i in range(len(u)):
        if i == 0:
            x[i] = a
        else:
            x[i] = x[i-1] + dx

        u_init[i] = np.math.sin(Pi*(x[i]+1))
        u_exact[i] = exact(1, x[i])

    # Plot exact solution
    fig, (ax1, ax2) = plt.subplots(2, sharex=True)
    ax1.plot(x,u_exact,label="Exact")
    ax1.set_ylabel("y")
    ax1.set_xlim([-1,1])
    ax2.set_xlabel("x")
    ax2.set_ylabel("error")

    # Set u
    for i in range(len(u)):
        u[i] = u_init[i]

    # Solve using FTCS
    while True:

        t += dt
        method = "FTCS"

        u = FTCS(dt, dx, u)

        if abs(t - 1) <= tol:
            ax1.plot(x,u,label=method)
            ax2.semilogy(x,calcError(t,x,u),label=method, color='g')
            break

    # Reset u
    for i in range(len(u)):
        u[i] = u_init[i]

    # Reset time
    t = 0.0

    # Solve using BTCS
    while True:

        t += dt
        method = "BTCS"

        u = BTCS(dt, dx, u)

        if abs(t - 1) <= tol:
            ax1.plot(x,u,label=method)
            ax2.semilogy(x,calcError(t,x,u),label=method, color='r')
            break

    # Reset u
    for i in range(len(u)):
        u[i] = u_init[i]

    # Reset time
    t = 0.0

    # Solve using CN
    while True:

        t += dt
        method = "CN"

        u = CN(dt, dx, u)

        if abs(t - 1) <= tol:
            ax1.plot(x,u,label=method)
            ax2.semilogy(x,calcError(t,x,u),label=method, color='c')
            break

    ax1.legend(loc=3)
    ax2.legend(loc=3)
    ax2.set_ylim([10**-18, 10**0])
    ax1.grid()
    ax2.grid()
    plt.savefig("Case%d" %(case) + ".pdf", bbox_inches='tight')
    plt.show()

##main(1)
##main(2)
##main(3)
##main(4)
##main(5)

a = np.zeros(4)
b = np.zeros(4)
c = np.zeros(4)
r = np.array([1.0,2.0,3.0,4.0])

a.fill(1)
b.fill(2)
c.fill(3)

print(a)
print(b)
print(c)
print(r)

print(tdma(a,b,c,r))