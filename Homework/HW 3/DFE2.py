import NumMethods as nm
import numpy as np

###########################################
""" Solve equations using direct methods """

def calcCoeffs(x):
    # Calculates the coefficients
    c = np.zeros(3)
    c[0] = -(x + 3.0)/(x + 1.0)     # a(x)
    c[1] = (x + 3.0)/(x + 1.0)**2.0   # b(x)
    c[2] = 2.0*(x + 1.0) + 3.0*c[1]   # f(x)

    return c

def direct():
    # Number of points
    N = 20

    # Boundary conditions
    y0a = 5.0
    y0b = 4.0

    # Domain bounaries
    a = 0.0
    b = 2.0

    # Length
    L = b - a

    # Stepsize
    h = L/(N+1)

    # Initialize tri-diagonal system
    a = np.zeros(N)
    b = np.zeros(N)
    c = np.zeros(N)
    r = np.zeros(N)

    # Left side bounary
    j = 0
    xj = h

    xRange = [0, xj]

    coeffs = calcCoeffs(xj)

    alpha = 1/(h**2.0) + coeffs[0]/(2.0*h)
    beta = coeffs[1] - 2.0/(h**2.0)
    gamma = 1/(h**2.0) - coeffs[0]/(2.0*h)

    a[j] = gamma
    b[j] = beta
    c[j] = alpha
    r[j] = coeffs[2] - gamma*y0a

    # Mid-range points j = 1 to N-1
    # range is end point non-inclusive, so this really stops at N-2
    for j in range(1, (N - 1)):

        xj = (j + 1)/(N+1)*L

        xRange.append(xj)

        coeffs = calcCoeffs(xj)

        alpha = 1/(h**2.0) + coeffs[0]/(2.0*h)
        beta = coeffs[1] - 2.0/(h**2.0)
        gamma = 1/(h**2.0) - coeffs[0]/(2.0*h)

        a[j] = gamma
        b[j] = beta
        c[j] = alpha
        r[j] = coeffs[2]

    # Last point
    j = N - 1
    xj = L - h

    xRange.append(xj)
    xRange.append(L)

    coeffs = calcCoeffs(xj)

    alpha = 1/(h**2.0) + coeffs[0]/(2.0*h)
    beta = coeffs[1] - 2.0/(h**2.0)
    gamma = 1/(h**2.0) - coeffs[0]/(2.0*h)

    a[j] = gamma
    b[j] = beta
    c[j] = alpha
    r[j] = coeffs[2] - alpha*y0b

    # Call Thomas algorithm to solve for temperatures
    T = np.empty([0])
    T = np.append(T,y0a)
    T = np.append(T, nm.tdma(a, b, c, r))
    T = np.append(T,y0b)

    return xRange, T
