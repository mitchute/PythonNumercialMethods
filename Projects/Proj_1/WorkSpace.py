

import NumMethods as nm
import numpy as np
import matplotlib.pyplot as plt


### Problem 1

def func(x):
    return np.cos(10*x)*np.sin(x)

def func_derivative(x):
    return (11*np.cos(11*x) - 9*np.cos(9*x)) / 2.0

def func_der_2(x):
    return (81*np.sin(9*x) - 121*np.sin(11*x)) / 2.0

#### Section 1

x0 = -1
xn = 1

x_p = np.linspace(x0, xn, 100)
f_x = func(x_p)

x_9 = np.array(np.linspace(x0, xn, 9))
f_9 = func(x_9)
p_9 = nm.lagrange(x_9, f_9, x_9)

x_17 = np.array(np.linspace(x0, xn, 17))
f_17 = func(x_17)
p_17 = nm.lagrange(x_17, f_17, x_17)

x_33 = np.array(np.linspace(x0, xn, 33))
f_33 = func(x_33)
p_33 = nm.lagrange(x_33, f_33, x_33)

x_65 = np.array(np.linspace(x0, xn, 65))
f_65 = func(x_65)
p_65 = nm.lagrange(x_65, f_65, x_65)


# f, ax1 = plt.subplots()

# ax1.plot(x_p, f_x, marker='', label='f(x)')
# ax1.plot(x_9, p_9, marker='x', linestyle='-', label="Lagrange")
# ax1.set_title('9 Points')
# ax1.set_ylim([-1,1])
# ax1.set_xlabel('x')
# ax1.set_ylabel('f(x)')
# lgd = plt.legend(loc='upper right')
# ax1.grid()

# plt.savefig("1-1-1.png", bbox_extra_artists=(lgd,), bbox_inches='tight')


# f, ax1 = plt.subplots()

# ax1.plot(x_p, f_x, marker='', label='f(x)')
# ax1.plot(x_17, p_17, marker='x', linestyle='-', label="Lagrange")
# ax1.set_title('17 Points')
# ax1.set_ylim([-1,1])
# ax1.set_xlabel('x')
# ax1.set_ylabel('f(x)')
# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# lgd = plt.legend(loc='upper right')
# ax1.grid()

# plt.savefig("1-1-2.png", bbox_extra_artists=(lgd,), bbox_inches='tight')


# f, ax1 = plt.subplots()

# ax1.plot(x_p, f_x, marker='', label='f(x)')
# ax1.plot(x_33, p_33, marker='x', linestyle='-', label="Lagrange")
# ax1.set_title('33 Points')
# ax1.set_ylim([-1,1])
# ax1.set_xlabel('x')
# ax1.set_ylabel('f(x)')
# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# lgd = plt.legend(loc='upper right')
# ax1.grid()

# plt.savefig("1-1-3.png", bbox_extra_artists=(lgd,), bbox_inches='tight')


# f, ax1 = plt.subplots()

# ax1.plot(x_p, f_x, marker='', label='f(x)')
# ax1.plot(x_65, p_65, marker='x', linestyle='-', label="Lagrange")
# ax1.set_title('65 Points')
# ax1.set_ylim([-1,1])
# ax1.set_xlabel('x')
# ax1.set_ylabel('f(x)')
# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# lgd = plt.legend(loc='upper right')
# ax1.grid()

# plt.savefig("1-1-4.png", bbox_extra_artists=(lgd,), bbox_inches='tight')

#### Section 2

# f_der = func_derivative(x_p)

# p_9_der = []
# p_17_der = []
# p_33_der = []
# p_65_der = []

# for i in x_9:
   # p_9_der.append(nm.lagrange_derivative(x_9, f_9, i))
# for i in x_17:
   # p_17_der.append(nm.lagrange_derivative(x_17, f_17, i))
# for i in x_33:
   # p_33_der.append(nm.lagrange_derivative(x_33, f_33, i))
# for i in x_65:
   # p_65_der.append(nm.lagrange_derivative(x_65, f_65, i))

# p_9_centralDer = nm.c2(x_9, f_9)
# p_17_centralDer = nm.c2(x_17, f_17)
# p_33_centralDer = nm.c2(x_33, f_33)
# p_65_centralDer = nm.c2(x_65, f_65)


# f, ax1 = plt.subplots()

# l1, = ax1.plot(x_p, f_der, label="f'(x)")
# l2, = ax1.plot(x_9, p_9_der, label="p'(x)", marker='x')
# l3, = ax1.plot(x_9, p_9_centralDer, label="CentralDiff", marker='^')
# ax1.set_title('9 Points')
# ax1.set_ylabel("Derivative")
# ax1.set_ylim([-10,10])
# lgd = ax1.legend((l1, l2, l3), ("f'(x)", "p'(x)", "Central Diff"), loc='upper right')

# plt.savefig("1-2-1.png", bbox_extra_artists=(lgd,), bbox_inches='tight')



# f, ax1 = plt.subplots()

# l1, = ax1.plot(x_p, f_der, label="f'(x)")
# l2, = ax1.plot(x_17, p_17_der, label="p'(x)", marker='x')
# l3, = ax1.plot(x_17, p_17_centralDer, label="CentralDiff", marker='^')
# ax1.set_title('17 Points')
# ax1.set_ylabel("Derivative")
# ax1.set_ylim([-10,10])
# lgd = ax1.legend((l1, l2, l3), ("f'(x)", "p'(x)", "Central Diff"), loc='upper right')

# plt.savefig("1-2-2.png", bbox_extra_artists=(lgd,), bbox_inches='tight')




# f, ax1 = plt.subplots()

# l1, = ax1.plot(x_p, f_der, label="f'(x)")
# l2, = ax1.plot(x_33, p_33_der, label="p'(x)", marker='x')
# l3, = ax1.plot(x_33, p_33_centralDer, label="CentralDiff", marker='^')
# ax1.set_title('33 Points')
# ax1.set_ylabel("Derivative")
# ax1.set_ylim([-10,10])
# lgd = ax1.legend((l1, l2, l3), ("f'(x)", "p'(x)", "Central Diff"), loc='upper right')

# plt.savefig("1-2-3.png", bbox_extra_artists=(lgd,), bbox_inches='tight')




# f, ax1 = plt.subplots()

# l1, = ax1.plot(x_p, f_der, label="f'(x)")
# l2, = ax1.plot(x_65, p_65_der, label="p'(x)", marker='x')
# l3, = ax1.plot(x_65, p_65_centralDer, label="CentralDiff", marker='^')
# ax1.set_title('65 Points')
# ax1.set_ylabel("Derivative")
# ax1.set_ylim([-10,10])
# lgd = ax1.legend((l1, l2, l3), ("f'(x)", "p'(x)", "Central Diff"), loc='upper right')

# plt.savefig("1-2-4.png", bbox_extra_artists=(lgd,), bbox_inches='tight')


### PROLBLEM 2

#### Section 1

# p_9_pade = nm.pade4(x_9, f_9)
# p_17_pade = nm.pade4(x_17, f_17)
# p_33_pade = nm.pade4(x_33, f_33)
# p_65_pade = nm.pade4(x_65, f_65)


# f, ax1 = plt.subplots()

# l1, = ax1.plot(x_p, f_der, label="f'(x)")
# l2, = ax1.plot(x_9, p_9_pade, label="Pade", marker='x')
# ax1.set_title('9 Points')
# ax1.set_ylabel("Derivative")
# ax1.set_ylim([-10,10])
# lgd = ax1.legend((l1, l2), ("f'(x)", "Pade"), loc='upper right')

# plt.savefig("2-1-1.png", bbox_extra_artists=(lgd,), bbox_inches='tight')


# f, ax1 = plt.subplots()

# l1, = ax1.plot(x_p, f_der, label="f'(x)")
# l2, = ax1.plot(x_17, p_17_pade, label="Pade", marker='x')
# ax1.set_title('17 Points')
# ax1.set_ylabel("Derivative")
# ax1.set_ylim([-10,10])
# lgd = ax1.legend((l1, l2), ("f'(x)", "Pade"), loc='upper right')

# plt.savefig("2-1-2.png", bbox_extra_artists=(lgd,), bbox_inches='tight')


# f, ax1 = plt.subplots()

# l1, = ax1.plot(x_p, f_der, label="f'(x)")
# l2, = ax1.plot(x_33, p_33_pade, label="Pade", marker='x')
# ax1.set_title('33 Points')
# ax1.set_ylabel("Derivative")
# ax1.set_ylim([-10,10])
# lgd = ax1.legend((l1, l2), ("f'(x)", "Pade"), loc='upper right')

# plt.savefig("2-1-3.png", bbox_extra_artists=(lgd,), bbox_inches='tight')


# f, ax1 = plt.subplots()

# l1, = ax1.plot(x_p, f_der, label="f'(x)")
# l2, = ax1.plot(x_65, p_65_pade, label="Pade", marker='x')
# ax1.set_title('65 Points')
# ax1.set_ylabel("Derivative")
# ax1.set_ylim([-10,10])
# lgd = ax1.legend((l1, l2), ("f'(x)", "Pade"), loc='upper right')

# plt.savefig("2-1-4.png", bbox_extra_artists=(lgd,), bbox_inches='tight')


####Section 2

def get_index(x):
   l = len(x)
   index = (l-1)*3/4
   return index

def percent_error(nm_soln, exact_soln):
   return abs((exact_soln - nm_soln)/exact_soln)*100

i_9 = int(get_index(x_9))
i_17 = int(get_index(x_17))
i_33 = int(get_index(x_33))
i_65 = int(get_index(x_65))

# pade_error = [percent_error(p_9_pade[i_9],func_derivative(0.5)),
             # percent_error(p_17_pade[i_17],func_derivative(0.5)),
             # percent_error(p_33_pade[i_33],func_derivative(0.5)),
             # percent_error(p_65_pade[i_65],func_derivative(0.5))]

# lagrange_error = [percent_error(p_9_der[i_9],func_derivative(0.5)),
                 # percent_error(p_17_der[i_17],func_derivative(0.5)),
                 # percent_error(p_33_der[i_33],func_derivative(0.5)),
                 # percent_error(p_65_der[i_65],func_derivative(0.5))]

# c2_error = [percent_error(p_9_centralDer[i_9],func_derivative(0.5)),
           # percent_error(p_17_centralDer[i_17],func_derivative(0.5)),
           # percent_error(p_33_centralDer[i_33],func_derivative(0.5)),
           # percent_error(p_65_centralDer[i_65],func_derivative(0.5))]

# h_error = [(x_9[1]-x_9[0]),
          # (x_17[1]-x_17[0]),
          # (x_33[1]-x_33[0]),
          # (x_65[1]-x_65[0])]

# f, ax1 = plt.subplots()

# l1, = ax1.loglog(h_error, pade_error, marker='o', label='Pade')
# l2, = ax1.loglog(h_error, lagrange_error, marker='x', label='Lagrange')
# l3, = ax1.loglog(h_error, c2_error, marker='x', label='C2')
# ax1.set_title('Error vs. step size')
# ax1.set_ylabel("% Error")
# ax1.set_xlabel("Step Size")
# ax1.grid(b=True, which='major', color='grey', linestyle='--')
# ax1.grid(b=True, which='minor', color='grey', linestyle='--')
# lgd = ax1.legend((l1, l2, l3), ("Pade", "Lagrange", "2nd-O Central"), loc='center right')

# plt.savefig("2-2-1.png", bbox_extra_artists=(lgd,), bbox_inches='tight')


### Problem 3

f_der_2 = func_der_2(x_p)

p_9_pade_2 = nm.pade4_2(x_9, f_9)
p_17_pade_2 = nm.pade4_2(x_17, f_17)
p_33_pade_2 = nm.pade4_2(x_33, f_33)
p_65_pade_2 = nm.pade4_2(x_65, f_65)

p_9_pade_using_pade4 = nm.pade4_twice(x_9, f_9)
p_17_pade_using_pade4 = nm.pade4_twice(x_17, f_17)
p_33_pade_using_pade4 = nm.pade4_twice(x_33, f_33)
p_65_pade_using_pade4 = nm.pade4_twice(x_65, f_65)


# f, ax1 = plt.subplots()

# l1, = ax1.plot(x_p, f_der_2, label='f"(x)')
# l2, = ax1.plot(x_9, p_9_pade_2, marker='x', label='Pade')
# ax1.set_title('9 Points')
# lgd = ax1.legend((l1, l2), ('f"(x)', "Pade"), loc='upper right')

# plt.savefig("3-1-1.png", bbox_extra_artists=(lgd,), bbox_inches='tight')


# f, ax1 = plt.subplots()

# l1, = ax1.plot(x_p, f_der_2, label='f"(x)')
# l2, = ax1.plot(x_17, p_17_pade_2, marker='x', label='Pade')
# ax1.set_title('17 Points')
# lgd = ax1.legend((l1, l2), ('f"(x)', "Pade"), loc='upper right')

# plt.savefig("3-1-2.png", bbox_extra_artists=(lgd,), bbox_inches='tight')



# f, ax1 = plt.subplots()

# l1, = ax1.plot(x_p, f_der_2, label='f"(x)')
# l2, = ax1.plot(x_33, p_33_pade_2, marker='x', label='Pade')
# ax1.set_title('33 Points')
# lgd = ax1.legend((l1, l2), ('f"(x)', "Pade"), loc='upper right')

# plt.savefig("3-1-3.png", bbox_extra_artists=(lgd,), bbox_inches='tight')



# f, ax1 = plt.subplots()

# l1, = ax1.plot(x_p, f_der_2, label='f"(x)')
# l2, = ax1.plot(x_65, p_65_pade_2, marker='x', label='Pade')
# ax1.set_title('65 Points')
# lgd = ax1.legend((l1, l2), ('f"(x)', "Pade"), loc='upper right')

# plt.savefig("3-1-4.png", bbox_extra_artists=(lgd,), bbox_inches='tight')



# f, ax1 = plt.subplots()

# l1, = ax1.plot(x_p, f_der_2, label='f"(x)')
# l2, = ax1.plot(x_9, p_9_pade_using_pade4, marker='x', label='Pade')
# ax1.set_title('9 Points')
# lgd = ax1.legend((l1, l2), ('f"(x)', "Pade"), loc='upper right')

# plt.savefig("3-2-1.png", bbox_extra_artists=(lgd,), bbox_inches='tight')



# f, ax1 = plt.subplots()

# l1, = ax1.plot(x_p, f_der_2, label='f"(x)')
# l2, = ax1.plot(x_17, p_17_pade_using_pade4, marker='x', label='Pade')
# ax1.set_title('17 Points')
# lgd = ax1.legend((l1, l2), ('f"(x)', "Pade"), loc='upper right')

# plt.savefig("3-2-2.png", bbox_extra_artists=(lgd,), bbox_inches='tight')


# f, ax1 = plt.subplots()

# l1, = ax1.plot(x_p, f_der_2, label='f"(x)')
# l2, = ax1.plot(x_33, p_33_pade_using_pade4, marker='x', label='Pade')
# ax1.set_title('33 Points')
# lgd = ax1.legend((l1, l2), ('f"(x)', "Pade"), loc='upper right')

# plt.savefig("3-2-3.png", bbox_extra_artists=(lgd,), bbox_inches='tight')



# f, ax1 = plt.subplots()

# l1, = ax1.plot(x_p, f_der_2, label='f"(x)')
# l2, = ax1.plot(x_65, p_65_pade_using_pade4, marker='x', label='Pade')
# ax1.set_title('65 Points')
# lgd = ax1.legend((l1, l2), ('f"(x)', "Pade"), loc='upper right')

# plt.savefig("3-2-4.png", bbox_extra_artists=(lgd,), bbox_inches='tight')




pade4_2_error = [percent_error(p_9_pade_2[i_9],func_der_2(0.5)),
             percent_error(p_17_pade_2[i_17],func_der_2(0.5)),
             percent_error(p_33_pade_2[i_33],func_der_2(0.5)),
             percent_error(p_65_pade_2[i_65],func_der_2(0.5))]

pade4_twice_error = [percent_error(p_9_pade_using_pade4[i_9],func_der_2(0.5)),
             percent_error(p_17_pade_using_pade4[i_17],func_der_2(0.5)),
             percent_error(p_33_pade_using_pade4[i_33],func_der_2(0.5)),
             percent_error(p_65_pade_using_pade4[i_65],func_der_2(0.5))]

h_error = [(x_9[1]-x_9[0]),
          (x_17[1]-x_17[0]),
          (x_33[1]-x_33[0]),
          (x_65[1]-x_65[0])]

f, ax1 = plt.subplots()

l1, = ax1.loglog(h_error, pade4_2_error, marker='o')
l2, = ax1.loglog(h_error, pade4_twice_error, marker='x')
ax1.set_title('Error vs. step size')
ax1.set_ylabel("% Error")
ax1.set_xlabel("Step Size")
ax1.grid(b=True, which='major', color='grey', linestyle='--')
ax1.grid(b=True, which='minor', color='grey', linestyle='--')
lgd = ax1.legend((l1, l2), ("Direct", "Indirect"), loc='center right')

plt.savefig("3-3-1.png", bbox_extra_artists=(lgd,), bbox_inches='tight')


### Problem 4

def abs_error(nm_soln, exact_soln):
    return abs(exact_soln - nm_soln)

# outFile = open("Integration_data.csv", 'w')

# trap_func = []
# simp_func = []
# trapEnd_func = []

# trap_func_error = []
# simp_func_error = []
# trapEnd_func_error = []

# for i in range(4):

   # if i == 0:
       # outFile.write("Num" + "," + "Trapezoidal" + "," + "Simpson's" + "," + "TrapEnd" + "\n")

   # n = [8, 16, 32, 64]
   # x_in = [x_9, x_17, x_33, x_65]
   # f_in = [f_9, f_17, f_33, f_65]

   # trap = nm.trap1D(x_in[i], f_in[i])
   # trap_func.append(abs(trap))

   # trap_error = abs_error(trap, 0.0)
   # trap_func_error.append(abs(trap_error))

   # simp = nm.simp1D(x_in[i], f_in[i])
   # simp_func.append(abs(simp))

   # simp_error = abs_error(simp, 0.0)
   # simp_func_error.append(abs(simp_error))

   # trapEnd = nm.trapEnd1D(x_in[i], f_in[i])
   # trapEnd_func.append(abs(trapEnd))

   # trapEnd_error = abs_error(trapEnd, 0.0)
   # trapEnd_func_error.append(abs(trapEnd_error))

   # outFile.write(str(n[i]) + "," + str(trap) + "," + str(simp) + "," + str(trapEnd) + "\n")
   # outFile.write(str(n[i]) + "," + str(trap_error) + "," + str(simp_error) + "," + str(trapEnd_error) + ",error\n")

# outFile.close()

# n = 4
# width = 0.25

# ind = np.arange(n)

# fig, ax1 = plt.subplots()

# rects1 = ax1.bar(ind, trap_func, width, color='r')
# rects2 = ax1.bar(ind+1*width, simp_func, width, color='y')
# rects3 = ax1.bar(ind+2*width, trapEnd_func, width, color='g')

# plt.gca().set_yscale('log')

# ax1.set_ylabel('Abs Error')
# ax1.set_xticks(ind+1.5*width)
# ax1.set_xticklabels( ('8', '16', '32', '64') )

# lgd = ax1.legend((rects1, rects2, rects3), ("Trap", "Simp", "TrapEnd"), loc='center right')

# plt.savefig("4-1.png", bbox_extra_artists=(lgd,), bbox_inches='tight')

x0 = 0
xn = 1

f_9_exp = []
f_17_exp = []
f_33_exp = []
f_65_exp = []

x_9_exp = np.linspace(x0, xn, 9)
x_17_exp = np.linspace(x0, xn, 17)
x_33_exp = np.linspace(x0, xn, 33)
x_65_exp = np.linspace(x0, xn, 65)

for i in x_9_exp:
    f_9_exp.append(np.math.exp(i))
for i in x_17_exp:
    f_17_exp.append(np.math.exp(i))
for i in x_33_exp:
    f_33_exp.append(np.math.exp(i))
for i in x_65_exp:
    f_65_exp.append(np.math.exp(i))

outFile = open("Integration_data_exp.csv", 'w')

trap_func = []
simp_func = []
trapEnd_func = []
gauss_func = []
exact_func = []

trap_func_error = []
simp_func_error = []
trapEnd_func_error = []
gauss_func_error = []

for i in range(4):

    if i == 0:
        outFile.write("Num" + "," + "Trapezoidal" + "," + "Simpson's" + "," + "TrapEnd" + "," + "Gauss" + "\n")

    n = [8, 16, 32, 64]
    x_in = [x_9_exp, x_17_exp, x_33_exp, x_65_exp]
    f_in = [f_9_exp, f_17_exp, f_33_exp, f_65_exp]

    trap = nm.trap1D(x_in[i], f_in[i])
    trap_func.append(trap)

    trap_error = abs_error(trap, 1.7182818284)
    trap_func_error.append(trap_error)

    simp = nm.simp1D(x_in[i], f_in[i])
    simp_func.append(simp)

    simp_error = abs_error(simp, 1.7182818284)
    simp_func_error.append(simp_error)

    trapEnd = nm.trapEnd1D(x_in[i], f_in[i])
    trapEnd_func.append(trapEnd)

    trapEnd_error = abs_error(trapEnd, 1.7182818284)
    trapEnd_func_error.append(trapEnd_error)

    gauss = nm.gaussQuadrature(0, 1)
    gauss_func.append(gauss)

    gauss_error = abs_error(gauss, 1.7182818284)
    gauss_func_error.append(gauss_error)

    exact_func.append(1.7182818284)

    outFile.write(str(n[i]) + "," + str(trap) + "," + str(simp) + ","
                  + str(trapEnd) + "," + str(gauss) + "\n")
    outFile.write(str(n[i]) + "," + str(trap_error) + "," + str(simp_error)
                  + "," + str(trapEnd_error) + "," + str(gauss_error)
                  + ",error\n")

outFile.close()

n = 4
width = 0.15

ind = np.arange(n)

fig, ax1 = plt.subplots()

rects1 = ax1.bar(ind, trap_func, width, color='r')
rects2 = ax1.bar(ind+1*width, simp_func, width, color='y')
rects3 = ax1.bar(ind+2*width, trapEnd_func, width, color='g')
rects4 = ax1.bar(ind+3*width, gauss_func, width, color='b')
rects5 = ax1.bar(ind+4*width, exact_func, width, color='orange')

ax1.set_ylabel('Integral')
ax1.set_xlabel('Num Data Points')
ax1.set_xticks(ind+2.5*width)
ax1.set_xticklabels( ('8', '16', '32', '64') )
ax1.set_ylim([1.55, 1.75])

lgd = ax1.legend((rects1, rects2, rects3, rects4, rects5), ("Trap", "Simp", "TrapEnd", "Gauss", "Exact"), loc='center right')

plt.savefig("4-2.png", bbox_extra_artists=(lgd,), bbox_inches='tight')

# fig, ax1 = plt.subplots()

# rects1 = ax1.bar(ind, trap_func_error, width, color='r')
# rects2 = ax1.bar(ind+1*width, simp_func_error, width, color='y')
# rects3 = ax1.bar(ind+2*width, trapEnd_func_error, width, color='g')
# rects4 = ax1.bar(ind+3*width, gauss_func_error, width, color='b')

# ax1.set_ylabel('Abs Error')
# ax1.set_xlabel('Num Data Points')
# ax1.set_xticks(ind+2*width)
# ax1.set_xticklabels( ('8', '16', '32', '64') )

# plt.gca().set_yscale('log')

# lgd = ax1.legend((rects1, rects2, rects3, rects4), ("Trap", "Simp", "TrapEnd", "Gauss"), loc='center right')

# plt.savefig("4-3.png", bbox_extra_artists=(lgd,), bbox_inches='tight')
