import matplotlib.pyplot as plt
import numpy as np
import NumMethods as nm

###########################################
""" Solve equations using shooting method """

# Number of equations
ne = 2

# Initializiations
y = np.zeros(ne)

# Number of  data points
n = 1000

# Initialize array for storing plot variables
y_plt = np.zeros((ne, n))

# Domain bounaries
a = 0.0
b = 2.0

# Stepsize
h = (b - a)/float(n)

# Boundary conditions
y0a = 5.0
y0b = 4.0

# Initial guesses for y[1] at x = 0
A1 = 1.0
A2 = 0.5

# First guess
# Initial condition
y[0] = y0a  # Boundary condition
y[1] = A1   # Guess

# Time integration with RK4 scheme to find B
for i in range(n):
    x = i*h
    y = nm.RK4(x,h,y)

# Assign estimated value for y[0] at x = 10
B1 = y[0]

# Second guess
# Initial condition
y[0] = y0a  # Boundary condition
y[1] = A2   # Guess

# Time integration with RK4 scheme
for i in range(n):
    x = i*h
    y = nm.RK4(x,h,y)

# Assign estimated value for y[0] at x = 10
B2 = y[0]

# Iteration for the rest, it will converge quickly
num_iter = 100
for i in range(num_iter):

    # Initial condition
    guess = A2 + (y0b - B2)/((B2 - B1)/(A2 - A1)) # Secant method

    y[0] = y0a      # Boundary condition
    y[1] = guess    # Guess

    # Time integration with RK4 scheme
    for j in range(n):
        x = j*h
        y = nm.RK4(x,h,y)

        # Update plot variables
        for k in range(ne):
            y_plt[k][j] = y[k]

    # Update for next iteration
    B1 = B2
    B2 = y[0]
    A1 = A2
    A2 = guess

    # Check for the final point
    tolerance = 10e-6
    if (np.abs(B2 - y0b) <= tolerance ):
        break

###########################################
""" Plot the solutions """

x_plt = np.arange(a,b,h)

fig, ax = plt.subplots()

ax.plot(x_plt, y_plt[0], label="T(x)")
ax.set_ylabel("Temp [C]")
ax.set_xlabel("x")

plt.tight_layout()
plt.grid()
plt.legend(loc='best')
plt.show()



