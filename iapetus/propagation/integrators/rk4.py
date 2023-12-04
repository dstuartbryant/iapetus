import math

import numpy as np


def rk4(funcptr, y0, tspan, dt, tspantol, extras=None):
    """Fixed-step 4th-order Runge-Kutta initial value ODE solver."""

    num_times = len(tspan)
    statesize = len(y0)
    Y = np.zeros((num_times, statesize))

    if tspan[0] < tspan[1]:
        dt = abs(dt)
    else:
        dt = -1.0 * abs(dt)

    # Determine number of subintervals of tspan
    itnum = np.zeros(num_times)
    for i in range(0, num_times - 1):
        frac = math.modf((tspan[i + 1] - tspan[i]) / dt)
        itnum[i] = frac[1]
        if math.fabs(frac[0]) > tspantol:
            itnum[i] += 1

    Y[0, :] = y0
    yi = y0
    y = np.zeros(statesize)
    for i in range(0, num_times - 1):
        tt = tspan[i]
        dtt = dt

        dt2 = dtt / 2.0
        dt6 = dtt / 6.0

        for j in range(0, int(itnum[i])):
            if (i == num_times - 2) and (j == itnum[i] - 1):
                dtt = tspan[i + 1] - tt

                dt2 = dtt / 2.0
                dt6 = dtt / 6.0

            k1 = funcptr(tt, yi, extras)
            for k in range(0, statesize):
                y[k] = yi[k] + dt2 * k1[k]

            k2 = funcptr(tt + dt2, y, extras)
            for k in range(0, statesize):
                y[k] = yi[k] + dt2 * k2[k]

            k3 = funcptr(tt + dt2, y, extras)
            for k in range(0, statesize):
                y[k] = yi[k] + dtt * k3[k]

            k4 = funcptr(tt + dtt, y, extras)
            for k in range(0, statesize):
                yi[k] = yi[k] + dt6 * (k1[k] + 2.0 * (k2[k] + k3[k]) + k4[k])

            tt = tt + dtt

        # end for j in ...

        Y[i + 1, :] = yi

    # end for i in ...

    return tspan, Y
