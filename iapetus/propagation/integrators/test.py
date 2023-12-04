import math

import numpy as np
from rk4 import rk4
from rk45 import rk45
from rk78 import rk78


def deriv_fcn(t, y0, extras=None):
    dy = [y0[1], 2 * math.tan(t) * (1 / math.cos(t)) ** 2]
    return dy


y0 = [0.0, 1.0]
dt = 0.05
tspan = np.arange(0, 1.4 + dt, dt)
tspantol = 2e-5


T4, Y4 = rk4(deriv_fcn, y0, tspan, dt, tspantol)

y0 = [0.0, 1.0]
T45, Y45 = rk45(deriv_fcn, y0, tspan, dt, tspantol)

y0 = [0.0, 1.0]
tspantol = 2e-10
T78, Y78 = rk78(deriv_fcn, y0, tspan, dt, tspantol)


print("\n\n")
print(
    "{:s}".format("time").ljust(7),
    "{:s}".format("rk4 result").ljust(20),
    "{:s}".format("rk45 result").ljust(20),
    "{:s}".format("rk78 result").ljust(20),
    "{:s}".format("tan(time)").ljust(20),
    "{:s}".format("rk4 error").ljust(20),
    "{:s}".format("rk45 error").ljust(20),
    "{:s}".format("rk78 error").ljust(20),
)
print("{:-^150}".format(""))
for idx in range(0, len(T4)):
    print(
        "{: 1.2f}".format(T4[idx]).ljust(7),
        "{: 1.11f}".format(Y4[idx][0]).ljust(20),
        "{: 1.11}".format(Y45[idx][0]).ljust(20),
        "{: 1.11}".format(Y78[idx][0]).ljust(20),
        "{: 1.11f}".format(math.tan(T4[idx])).ljust(20),
        "{: 1.11f}".format(math.tan(T4[idx]) - Y4[idx][0]).ljust(20),
        "{: 1.11f}".format(math.tan(T45[idx]) - Y45[idx][0]).ljust(20),
        "{: 1.11f}".format(math.tan(T78[idx]) - Y78[idx][0]).ljust(20),
    )
