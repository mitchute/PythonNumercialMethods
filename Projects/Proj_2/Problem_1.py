import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

################################################################################

sigma = 10.0
beta = 8.0/3.0
rho = 28.0

################################################################################

def Euler(x,h,y):

    y_out = h * RHS(x,y) + y

    return y_out

################################################################################

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

################################################################################

def RHS(x,y):
    # Right hand side
    n = len(y)
    r = np.zeros(n)

    r[0] = sigma * (y[1] - y[0])
    r[1] = y[0] * (rho - y[2]) - y[1]
    r[2] = y[0] * y[1] - beta * y[2]

    return r

################################################################################

def main(method,mesh):
    # Number of equations
    ne = 3

    # Domain bounaries
    a = 0.0
    b = 100.0

    # Stepsize
    if mesh == "coarse":
        h = 0.01
        h_text = "coarse"
    elif mesh == "fine":
        h = 0.001
        h_text = "fine"

    # Initial conditions
    y0a = 1.0
    y1a = 1.0
    y2a = 1.0

    # Initialize
    y = np.array([y0a, y1a, y2a])
    time = np.arange(a,b+h,h)
    colors = np.linspace(0,1,len(time))

    # Initialize array for storing plot variables
    y_plt = np.zeros((ne, len(time)))

    # Time integration with RK4 scheme
    for i in range(len(time)):
        x = 1 #Dummy variable

        if method == "Euler":
            y = Euler(x,h,y)
        elif method == "RK":
            y = RK4(x,h,y)

        # Update plot variables
        for j in range(ne):
            y_plt[j][i] = y[j]

    # Make Plots

    # x-t Plot
    fig, ax = plt.subplots()
    ax.plot(time,y_plt[0], color='r')
    ax.set_xlabel("t")
    ax.set_ylabel("x")
    plt.grid()
    plt.savefig("x-" + method + "-" + h_text +".pdf", bbox_inches='tight')

    # y-t Plot
    fig, ax = plt.subplots()
    ax.plot(time,y_plt[1], color='g')
    ax.set_xlabel("t")
    ax.set_ylabel("y")
    plt.grid()
    plt.savefig("y-" + method + "-" + h_text +".pdf", bbox_inches='tight')

    # z-t Plot
    fig, ax = plt.subplots()
    ax.plot(time,y_plt[2], color='b', )
    ax.set_xlabel("t")
    ax.set_ylabel("z")
    plt.grid()
    plt.savefig("z-" + method + "-" + h_text +".pdf", bbox_inches='tight')

    # 3-D Plot
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    for i in range(len(y_plt[0])):
        ax.plot(y_plt[0,i:i+2], y_plt[1,i:i+2], y_plt[2,i:i+2], color=plt.cm.jet(int(255*i/len(y_plt[0]))))
    plt.savefig("3D-" + method + "-" + h_text +".pdf", bbox_inches='tight')

main("Euler", "coarse")
main("Euler", "fine")
main("RK", "coarse")
main("RK", "fine")

