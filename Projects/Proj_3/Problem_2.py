import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

################################################################################

# Gauss-Seidel Solution (GS)
def GS(nx,ny,dx,dy,tol):

    numIter = 0
    error = 100 # initalized big, will get updated

    v = np.zeros((nx+1,ny+1))
    u = np.zeros((nx+1,ny+1))

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

                    v[i,j] = np.math.sin(2*np.math.pi*y)

                elif (j == 0 or j == ny):
                # Lower and upper boundary

                    v[i,j] = 0

                elif (i == nx):
                # Right-side boundary

                    v[i,j] = 0

                else:
                # Central points
                    f = 30*(x*x - x) + 30*(y*y - y)
                    v[i,j] = (1./4.)*(u[i+1,j] + u[i-1,j] + u[i,j+1] + u[i,j-1]
                    - (dx**2)*f)

        # Calculate error
        error = calcError(u,v)

        # Increment iteration counter
        numIter += 1

        # Save array for error calculation
        for j in range(ny+1):
            for i in range(nx+1):
                u[i,j] = v[i,j]

    return u, numIter, "Gauss-Seidel", "GS"

################################################################################

def calcError(arr1,arr2):

    e = np.zeros((len(arr1),len(arr2)))

    # Calculate error
    e = np.abs(np.subtract(arr1,arr2))

    # Biggest error
    error = np.amax(e)

    return error

################################################################################

def makePlot(x,y,z,name,error):
    # Setup x and y data ranges for plotting
    x_plt, y_plt = np.meshgrid(x, y, indexing='ij')

    # Plot solution
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_surface(x_plt, y_plt, z, cmap=plt.cm.jet, rstride=1, cstride=1, linewidth=0)

    # Equal aspect ratio
    ax.set_aspect('equal')

    # Set plot orientation
    ax.view_init(azim=-125, elev=25)

    # Plot formatting
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel(r'$\phi$')

    if name != "Exact":
        label = "Max Error: %0.3e" %(error)
        ax.text(0.1, 0.5, -1, label, None)

    # Set axis limit
    ax.set_zlim([-1,1])

    # Figure name
    figName = "%s.pdf" %(name)

    # Save/show figure
    plt.savefig(figName, bbox_inches='tight')
    plt.show()

################################################################################

def exact(x,y):

    z = np.zeros((len(x),len(y)))

    PI = np.math.pi

    for j in range(len(y)):
        for i in range(len(x)):
            term1 = 15*(x[i]**2 - x[i])*(y[j]**2 - y[j])
            term2 = np.math.sin(2*PI*y[j])
            term3 = np.math.sinh(2*PI*(x[i]-1))
            term4 = np.math.sinh(2*PI)
            z[i,j] = term1 - term2 * term3 / term4
    return z

################################################################################

# Discrete Sine Transform
def DST(nx,ny,dx,dy):

    N = nx+1
    M = ny+1

    gt = np.zeros((N,M))
    ft = np.zeros((N,M))
    f = np.zeros((N,M))
    u = np.zeros((N,M))

    PI = np.math.pi

    # 1: Fast inverse discrete sine transform of source term
    for k in range(1,nx):
        for l in range(1,ny):

            # sum
            for i in range(1,nx):
                x = i*dx

                for j in range(1,ny):
                    y = j*dy

                    term1 = 30*(x*x - x) + 30*(y*y - y)
                    term2 = (4*PI**2)*(x - 1)*np.math.sin(2*PI*y)

                    g = term1 - term2

                    term3 = np.math.sin(PI*float(k)*float(i)/float(nx))
                    term4 = np.math.sin(PI*float(l)*float(j)/float(ny))

                    gt[k,l] += g*term3*term4

    # Normalize
    for k in range(1,nx):
        for l in range(1,ny):
            gt[k,l] = gt[k,l]*(2./float(N))*(2./float(M))

    # 2: Compute Fourier coefficient of solution f:
    for k in range(1,nx):
        for l in range(1,ny):
            term1 = (2./(dx*dx))*(np.math.cos(PI*float(k)/float(nx)) - 1.)
            term2 = (2./(dy*dy))*(np.math.cos(PI*float(l)/float(ny)) - 1.)
            alpha = term1 + term2

            ft[k,l] = gt[k,l]/alpha

    # 3: Fast forward Fourier sine transform to find f:
    for i in range(1,nx):
        for j in range(1,ny):

            # sum
            for k in range(1,nx):
                for l in range(1,ny):
                    term1 = np.math.sin(PI*float(k)*float(i)/float(nx))
                    term2 = np.math.sin(PI*float(l)*float(j)/float(ny))

                    f[i,j] += ft[k,l]*term1*term2

    # 4: Transform by removing introduced varialbe to find u:
    for i in range(N):
        x = i*dx

        for j in range(M):
            y = j*dy
            u[i,j] = f[i,j] - (x - 1)*np.math.sin(2*PI*y)

    return u

################################################################################

##### Common variables #####

x0 = 0
xL = 1

y0 = 0
yL = 1

# Tolerance
tol = 10**(-10)

##### Numerical Solution #####

# Number of points
nx = 31
ny = 31

# Stepsize
dx = (xL - x0)/float(nx)
dy = (yL - y0)/float(ny)

# setup x, y, and u ranges
x_in = np.arange(x0,xL+dx,dx)
y_in = np.arange(y0,yL+dy,dy)

# Get exact solution
u_e = exact(x_in, y_in)

# Make plot for exact solution
makePlot(x_in, y_in, u_e,"Exact", 0.0)

# Call numerical method
u_nm, numIter, methodName, methodAbrev = GS(nx,ny,dx,dy,tol)

# Print GS statistics
print("Num. Iteration for GS: %d" %(numIter))
print("Max Error for GS: %0.3e" %(calcError(u_e,u_nm)))

# Make plot for numerical method
makePlot(x_in, y_in, u_nm, "GS", calcError(u_e,u_nm))

# Call DST
u_dst = DST(nx,ny,dx,dy)

# Print DST statistics
print("Max Error for DST: %0.3e" %(calcError(u_e,u_dst)))

# Make plot for DST solution
makePlot(x_in, y_in, u_dst, "DST", calcError(u_e,u_dst))
