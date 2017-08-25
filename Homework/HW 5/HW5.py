import matplotlib.pyplot as plt
import numpy as np

# Exact solution
def u(x,y):
    numIter = 100
    a = 0.0
    PI = np.math.pi

    for n in range(1, (numIter), 2):

        frac = 1 / ((n*PI)**2 * np.math.sinh(2*n*PI))
        term1 = np.math.sinh(n*PI*x)
        term2 = np.math.cos(n*PI*y)
        a += frac * term1 * term2

    u = x/4.0 - 4.0*a

    return u

# Calculate error
def calcError(xRange,yRange,uRange):

    maxErr = 0.0
    percentErr = 0.0

    for j in range(len(yRange)):
        y = yRange[j]
        for i in range(len(xRange)):
            # Skip left and right boundary points
            if i == 0 or i == len(xRange)-1:
                continue
            x = xRange[i]

            analytical = u(x,y)
            numerical = uRange[i,j]

            e = abs(analytical - numerical)

            if (analytical != 0):
                p = abs(e/analytical)
                if p > percentErr:
                    percentErr = p

            if e > maxErr:
                maxErr = e

    return maxErr, percentErr

# Point Jacobi (PJ)
def PJ(nx,ny,dx,dy,tol):

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

                # Left side boundary
                if i == 0:
                    u[i,j] = 0
                    continue

                # Right side boundary
                elif i == nx:
                    u[i,j] = y
                    continue

                frac = 1/((-2.0/dx**2)+(-2.0/dy**2))

                # Lower boundary
                if j == 0:
                    u[i,j] = frac * ((-v[i+1,j] - v[i-1,j])/dx**2 - 2*v[i,j+1]/dy**2)
                    continue

                # Upper boundary
                elif j == ny:
                    u[i,j] = frac * ((-v[i+1,j] - v[i-1,j])/dx**2 - 2*v[i,j-1]/dy**2)
                    continue

                # Central points
                u[i,j] = frac * ((-v[i+1,j] - v[i-1,j])/dx**2 + (-v[i,j+1] - v[i,j-1])/dy**2)

        # Calculate error
        e = np.abs(np.subtract(u,v))

        # Biggest error
        error = np.amax(e)
        numIter += 1

        # Save array
        for j in range(ny+1):
            for i in range(nx+1):
                v[i,j] = u[i,j]

    return u, numIter, "Point-Jacobi", "PJ"

# Gauss-Seidel (GS)
def GS(nx,ny,dx,dy,tol):

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

                # Left side boundary
                if i == 0:
                    u[i,j] = 0
                    continue

                # Right side boundary
                elif i == nx:
                    u[i,j] = y
                    continue

                frac = 1/((-2.0/dx**2)+(-2.0/dy**2))

                # Lower boundary
                if j == 0:
                    u[i,j] = frac * ((-u[i+1,j] - u[i-1,j])/dx**2 - 2*u[i,j+1]/dy**2)
                    continue

                # Upper boundary
                elif j == ny:
                    u[i,j] = frac * ((-u[i+1,j] - u[i-1,j])/dx**2 - 2*u[i,j-1]/dy**2)
                    continue

                # Central points
                u[i,j] = frac * ((-u[i+1,j] - u[i-1,j])/dx**2 + (-u[i,j+1] - u[i,j-1])/dy**2)

        # Calculate error
        e = np.abs(np.subtract(u,v))

        # Biggest error
        error = np.amax(e)
        numIter += 1

        # Save array
        for j in range(ny+1):
            for i in range(nx+1):
                v[i,j] = u[i,j]

    return u, numIter, "Gauss-Seidel", "GS"

# Sucessive over relaxation (SOR)
def SOR(nx,ny,dx,dy,tol):

    numIter = 0
    error = 100 # initalized big, will get updated
    omega = 1.8

    up = 0.0
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

                # Left side boundary
                if i == 0:
                    up = 0
                    u[i,j] = u[i,j] + omega*(up-u[i,j])
                    continue

                # Right side boundary
                elif i == nx:
                    up = y
                    u[i,j] = u[i,j] + omega*(up-u[i,j])
                    continue

                frac = 1/((-2.0/dx**2)+(-2.0/dy**2))

                # Lower boundary
                if j == 0:
                    up = frac * ((-u[i+1,j] - u[i-1,j])/dx**2 - 2*u[i,j+1]/dy**2)
                    u[i,j] = u[i,j] + omega*(up-u[i,j])
                    continue

                # Upper boundary
                elif j == ny:
                    up = frac * ((-u[i+1,j] - u[i-1,j])/dx**2 - 2*u[i,j-1]/dy**2)
                    u[i,j] = u[i,j] + omega*(up-u[i,j])
                    continue

                # Central points
                up = frac * ((-u[i+1,j] - u[i-1,j])/dx**2 + (-u[i,j+1] - u[i,j-1])/dy**2)
                u[i,j] = u[i,j] + omega*(up-u[i,j])

        # Calculate error
        e = np.abs(np.subtract(u,v))

        # Biggest error
        error = np.amax(e)
        numIter += 1

        # Save array
        for j in range(ny+1):
            for i in range(nx+1):
                v[i,j] = u[i,j]

    return u, numIter, "Successive Over Relaxation", "SOR"

##### Common variables #####
x0 = 0
xL = 2

y0 = 0
yL = 1

tol = 0.00005

##### Setup plot for exact solution #####
step = 0.01

x_exact = np.arange(x0,xL+step,step)
y_exact = np.arange(y0,yL+step,step)
u_exact = np.zeros((len(x_exact),len(y_exact)))

# Fill the exact solution u-array
for i in range(x_exact.shape[0]):
    x = x_exact[i]
    for j in range(y_exact.shape[0]):
        y = y_exact[j]
        u_exact[i,j] = u(x,y)

# Setup x and y data ranges for plotting
x_exact, y_exact = np.meshgrid(x_exact, y_exact, indexing='ij')

##### Numerical Solution #####

# Number of points
nx = 10
ny = 10

# Stepsize
dx = (xL - x0)/float(nx)
dy = (yL - y0)/float(ny)

# setup x, y, and u ranges
x_nm = np.arange(x0,xL+dx,dx)
y_nm = np.arange(y0,yL+dy,dy)

# Call numerical method
#u_nm, numIter, methodName, methodAbrev = PJ(nx,ny,dx,dy,tol)
#u_nm, numIter, methodName, methodAbrev = GS(nx,ny,dx,dy,tol)
u_nm, numIter, methodName, methodAbrev = SOR(nx,ny,dx,dy,tol)

# Calculate error
maxError, percentError = calcError(x_nm, y_nm,u_nm)

# Setup x and y data ranges for plotting
x_nm, y_nm = np.meshgrid(x_nm, y_nm, indexing='ij')

# Plot solution
fig, ax = plt.subplots()
ax.contour(x_nm, y_nm, u_nm, 16, linestyles='dashed', colors='k', linewidths=1.5)
CS = ax.contour(x_exact, y_exact, u_exact, 16, linewidths=1.5)

# Equal aspect ratio
ax.set_aspect('equal')

# Set grid lines
horizTicks = np.arange(x0,xL+dx,dx)
verTicks = np.arange(y0,yL+dy,dy)
ax.set_xticks(horizTicks, minor=True)
ax.set_yticks(verTicks, minor=True)

ax.grid(which='minor')

# Add statistics to plot
ax.text(0.05,0.05,"Num Iter: %d" %(numIter),
        bbox={'facecolor':'white', 'alpha':1, 'pad':10, 'edgecolor':'black'})
ax.text(0.55,0.05,"Max Error: %0.2e" %(maxError) + r'$^\circ$' + 'C',
        bbox={'facecolor':'white', 'alpha':1, 'pad':10, 'edgecolor':'black'})
ax.text(1.35,0.05,"Per Error: %0.2f" %(percentError*100) + "%",
        bbox={'facecolor':'white', 'alpha':1, 'pad':10, 'edgecolor':'black'})

# Add color bar
cbar = plt.colorbar(CS, fraction=0.0234, pad=0.04)
cbar.ax.set_ylabel('Temperature')

# Plot formatting
plt.title("%s" %(methodName) + "\n" + r'$N_x = N_y = %d$' %(nx+1))
plt.ylabel('y')
plt.xlabel('x')

# Figure name
figName = "%s_%d.pdf" %(methodAbrev, nx+1)

# Save/show figure
plt.savefig(figName, bbox_inches='tight')
plt.show()

# Print stats
print("Num Iter: " + str(numIter))
print("Max Error: " + str(maxError))
print("Abs Error: " + str(percentError))
