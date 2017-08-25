import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Sucessive over relaxation (SOR)
def SOR(nx,ny,dx,dy,tol,omega):

    numIter = 0
    error = 100 # initalized big, will get updated

    v = np.zeros((nx+1,ny+1))
    u = np.zeros((nx+1,ny+1))
    e = np.zeros((nx+1,ny+1))

    while error >= tol:

        # Update u
        for j in range(ny+1):
            # y location
            y = j*dy

            for i in range(nx+1):
                # x location
                x = i*dx

                if (i == 0):
                # Left-side boundary

                    u[i,j] = y

                elif (j == 0):
                # Lower boundary

                    u[i,j] = 0

                elif ((i >= (nx/2) and j >= (ny/2)) or (j == ny)):
                # Middle step and upper boundary

                    u[i,j] = 1

                elif (i == nx):
                # Right-side boundary

                    u[i,j] = (1./11.)*(18.*u[-2,j] - 9*u[-3,j] +2*u[-4,j])

                else:
                # Central points

                    u[i,j] = u[i,j] + (omega/4.)*(u[i+1,j] + u[i-1,j] + u[i,j+1] + u[i,j-1] - 4*u[i,j])

        # Calculate error
        e = np.abs(np.subtract(u,v))

        # Biggest error
        error = np.amax(e)
        numIter += 1

        # Save array for error calculation
        for j in range(ny+1):
            for i in range(nx+1):
                v[i,j] = u[i,j]

    return u, numIter, "Successive Over Relaxation", "SOR"

def main(omega):

    ##### Common variables #####
    x0 = 0
    xL = 2

    y0 = 0
    yL = 1

    # Tolerance
    tol = 10**(-10)

    ##### Numerical Solution #####

    # Number of points
    nx = 20
    ny = 10

    # Stepsize
    dx = (xL - x0)/float(nx)
    dy = (yL - y0)/float(ny)

    # setup x, y, and u ranges
    x_nm = np.arange(x0,xL+dx,dx)
    y_nm = np.arange(y0,yL+dy,dy)

    # Call numerical method
    u_nm, numIter, methodName, methodAbrev = SOR(nx,ny,dx,dy,tol,omega)

    # Setup x and y data ranges for plotting
    x_nm, y_nm = np.meshgrid(x_nm, y_nm, indexing='ij')

    # Plot solution
    fig, ax = plt.subplots()
    CS = ax.contour(x_nm, y_nm, u_nm, 15, linewidths=1.5)

    # Equal aspect ratio
    ax.set_aspect('equal')

    # Set grid lines
    horizTicks = np.arange(x0,xL+dx,dx)
    verTicks = np.arange(y0,yL+dy,dy)
    ax.set_xticks(horizTicks, minor=True)
    ax.set_yticks(verTicks, minor=True)

    ax.grid(which='minor')

    # Add color bar
    cbar = plt.colorbar(CS, fraction=0.0234, pad=0.04)
    cbar.ax.set_ylabel('stream-function')

    # Plot formatting
    plt.ylabel('y')
    plt.xlabel('x')

    # Figure name
    figName = "omega-%d.pdf" %(omega*100)

    # Add rectangle
    ax.add_patch(mpatches.Rectangle((1.0,0.5), 1.0, 0.5, facecolor='grey', alpha=1.0, linewidth=2))

    # Add statistics
    props = dict(boxstyle='round', facecolor='white')
    ax.text(0.05, 0.95, " Num. Iter: %d \n Omega: %0.2f" %(numIter,omega), transform=ax.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)

    # Save/show figure
    plt.savefig(figName, bbox_inches='tight')

    return numIter

# Plot omega vs. iters
stepSize = 0.05
omegas = []
iters = []

for i in np.arange(1.0,2.0,stepSize):
    omegas.append(i)
    iters.append(main(i))

fig, ax = plt.subplots()
ax.plot(omegas, iters)
ax.grid()

minValIndex = iters.index(min(iters))

# Add statistics
props = dict(boxstyle='round', facecolor='white')
ax.text(0.05, 0.95, " Min Iter: %d \n Omega: %0.2f" %(iters[minValIndex],omegas[minValIndex]), transform=ax.transAxes, fontsize=14,
    verticalalignment='top', bbox=props)

# Plot formatting
ax.set_xlabel(r'$\omega$')
ax.set_ylabel("Iterations")
ax.set_ylim([0,300])
plt.savefig("NumIter.pdf", bbox_inches='tight')

plt.show()