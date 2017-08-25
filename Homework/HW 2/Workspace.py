import matplotlib.pyplot as plt
import numpy as np
import NumMethods as nm

# Number of equations
ne = 3

# Initializiations
y = np.zeros(ne)

# Number of  data points
n = 1000

# Initialize array for storing plot variables
y_plt = np.zeros((ne, n))

# Domain bounaries
a = 0.0
b = 10.0

# Stepsize
h = (b - a)/float(n)

# Boundary conditions
y2a = 0.0
y3a = 0.0
y2b = 1.0

# Initial guesses for y[1] at x = 0
A1 = 1.0
A2 = 0.5

# First guess
# Initial condition
y[0] = A1   # Guess
y[1] = y2a  # Boundary condition
y[2] = y3a  # Boundary condition

# Time integration with RK4 scheme to find B
for i in range(n):
    x = i*h
    y = nm.RK4(x,h,y)

# Assign estimated value for y[1] at x = 10
B1 = y[1]

# Second guess
# Initial condition
y[0] = A2   # Guess
y[1] = y2a  # Boundary condition
y[2] = y3a  # Boundary condition

# Time integration with RK4 scheme
for i in range(n):
    x = i*h
    y = nm.RK4(x,h,y)

# Assign estimated value for y[1] at x = 10
B2 = y[1]

# Iteration for the rest, it will converge quickly
num_iter = 100
for i in range(num_iter):

    # Initial condition
    guess = A2 + (y2b - B2)/((B2 - B1)/(A2 - A1)) # Secant method

    y[0] = guess    # Guess
    y[1] = y2a      # Boundary condition
    y[2] = y3a      # Boundary condition

    # Time integration with RK4 scheme
    for j in range(n):
        x = j*h
        y = nm.RK4(x,h,y)

        # Update plot variables
        for k in range(ne):
            y_plt[k][j] = y[k]

    # Update for next iteration
    B1 = B2
    B2 = y[1]
    A1 = A2
    A2 = guess

    # Check for the final point
    tolerance = 10e-6
    if (np.abs(B2 - y2b) <= tolerance ):
        break

# Plot
x_plt = np.arange(a,b,h)
x_ticks = np.arange(0,8,2)
y_ticks = np.arange(0,2.25,0.25)

fig, ax = plt.subplots()

ax.plot(x_plt, y_plt[0], label="f ''")
ax.plot(x_plt, y_plt[1], label="f '", linestyle='--')
ax.plot(x_plt, y_plt[2], label="f", linestyle=":")
ax.set_xlim([0,6])
ax.set_ylim([0,2])
ax.set_xticks(x_ticks)
ax.set_yticks(y_ticks)
ax.set_xlabel(r'$\eta$', fontsize=20)
ax.set_aspect('equal')

plt.tight_layout()
plt.grid()
plt.legend(loc='best')
plt.savefig("Blasius.pdf",bbox_inches='tight')



