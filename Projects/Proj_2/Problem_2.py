import numpy as np
import matplotlib.pyplot as plt

################################################################################

def exact(x,t):
    return np.math.sin(2*np.math.pi*(x-t))

################################################################################

def calcError(x,t,u):

    r = np.zeros(len(u))

    for i in range(len(u)):
        ue = exact(x[i],t)
        r[i] =  abs(ue - u[i])

    return r
################################################################################

def Euler(c,u):

    r = np.zeros(len(u))

    for i in range(len(u)):

        if i == 0:
            r[i] = u[i] - c * (u[i] - u[-2])
        else:
            r[i] = u[i] - c * (u[i] - u[i-1])

    return r

################################################################################

def Lax(c,u):

    r = np.zeros(len(u))

    for i in range(len(u)):

        term1 = 0.0
        term2 = 0.0

        if i == 0:

            term1 = (u[i+1] - u[-2])
            term2 = (u[i+1] - 2*u[i] + u[-2])
            r[i] = u[i] - (c/2.0) * term1 + ((c**2.0)/2.0) * term2

        elif i == len(u) - 1:

            term1 = (u[1] - u[i-1])
            term2 = (u[1] - 2*u[i] + u[i-1])
            r[i] = u[i] - (c/2.0) * term1 + ((c**2.0)/2.0) * term2

        else:

            term1 = (u[i+1] - u[i-1])
            term2 = (u[i+1] - 2*u[i] + u[i-1])
            r[i] = u[i] - (c/2.0) * term1 + ((c**2.0)/2.0) * term2

    return r

################################################################################

def main(c, method):

    a = 0.0
    b = 1.0
    h = 0.025
    dt = c * h
    t = 0.0
    tol = 0.00005

    # Initializations
    x = np.arange(a,b+h,h)
    u = np.zeros(len(x))

    for i in range(len(x)):
        u[i] = np.math.sin(2*np.math.pi*x[i])

    # Plot initial conditions
    fig, ax = plt.subplots()
    ax.plot(x,u, label="t=0")
    ax.set_xlabel("t")
    ax.set_ylabel("y")

    while True:

        t += dt

        if method == "Euler":
            u = Euler(c,u)
        elif method == "Lax":
            u = Lax(c,u)

        if abs(t - 1) <= tol:
            ax.plot(x,u,label="t=1")
            x1error = np.array(x)
            u1error = np.array(calcError(x,t,u))

        elif abs(t - 10) <= tol:
            ax.plot(x,u,label="t=10")
            x10error = np.array(x)
            u10error = np.array(calcError(x,t,u))
            break

    # Final plot
    plt.legend()
    plt.grid()
    plt.savefig(method + "-" + str(int(c*10)) + ".pdf", bbox_inches='tight')
    plt.show()


    # Error plot
    fig, ax = plt.subplots()
    ax.set_xlabel("t")
    ax.set_ylabel("error")
    ax.semilogy(x1error,u1error,label="t=1",color='g')
    ax.semilogy(x10error,u10error,label="t=10",color='r')
    ax.set_ylim([10**-18, 10**0])
    plt.legend()
    plt.grid()
    plt.savefig(method + "-error-" + str(int(c*10)) + ".pdf", bbox_inches='tight')

main(0.5, "Euler")
main(1.0, "Euler")
main(0.5, "Lax")
main(1.0, "Lax")
