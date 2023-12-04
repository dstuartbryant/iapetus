import math

import numpy as np

EPS = np.spacing(1)


def rk45(funcptr, y0, tspan, dt, tspantol, extras=None):
    dt0 = dt
    tol = tspantol

    namestr = "rk45"

    statesize = len(y0)
    tspansize = len(tspan)
    Y = np.zeros((tspansize, statesize))
    T = np.zeros((tspansize,))
    f = np.zeros((statesize * 7,))
    X = np.zeros((statesize,))
    dX = np.zeros((statesize,))
    X4 = np.zeros((statesize,))
    X5 = np.zeros((statesize,))
    Xcurrent = np.zeros((statesize,))
    Xcurrent5 = np.zeros((statesize,))

    # Set coefficient matrices
    alpha = np.zeros((6,))
    alpha[0] = 1.0 / 5.0
    alpha[1] = 3.0 / 10.0
    alpha[2] = 4.0 / 5.0
    alpha[3] = 8.0 / 9.0
    alpha[4] = 1.0
    alpha[5] = 1.0

    beta = np.zeros((6, 6))
    beta[0][0] = 1.0 / 5.0
    beta[0][1] = 3.0 / 40.0
    beta[0][2] = 44.0 / 45.0
    beta[0][3] = 19372.0 / 6561.0
    beta[0][4] = 9017.0 / 3168.0
    beta[0][5] = 35.0 / 384.0
    beta[1][1] = 9.0 / 40.0
    beta[1][2] = -56.0 / 15.0
    beta[1][3] = -25360.0 / 2187.0
    beta[1][4] = -355.0 / 33.0
    beta[2][2] = 32.0 / 9.0
    beta[2][3] = 64448.0 / 6561.0
    beta[2][4] = 46732.0 / 5247.0
    beta[2][5] = 500.0 / 1113.0
    beta[3][3] = -212.0 / 729.0
    beta[3][4] = 49.0 / 176.0
    beta[3][5] = 125.0 / 192.0
    beta[4][4] = -5103.0 / 18656.0
    beta[4][5] = -2187.0 / 6784.0
    beta[5][5] = 11.0 / 84.0

    chi4 = np.zeros((7,))
    chi4[0] = 5179.0 / 57600.0
    chi4[1] = 0.0
    chi4[2] = 7571.0 / 16695.0
    chi4[3] = 393.0 / 640.0
    chi4[4] = -92097.0 / 339200.0
    chi4[5] = 187.0 / 2100.0
    chi4[6] = 1.0 / 40.0

    chi5 = np.zeros((7,))
    chi5[0] = 35.0 / 384.0
    chi5[1] = 0.0
    chi5[2] = 500.0 / 1113.0
    chi5[3] = 125.0 / 192.0
    chi5[4] = -2187.0 / 6784.0
    chi5[5] = 11.0 / 84.0
    chi5[6] = 0.0

    # Set initial points in the output variables
    for i in range(0, statesize):
        Y[0][i] = y0[i]
    T[0] = tspan[0]

    # Initialize f
    for i in range(0, statesize):
        for j in range(0, 7):
            f[i + j * statesize] = 0.0

    time = tspan[0]
    htol = abs(time) * EPS
    sign = 1.0
    if tspan[tspansize - 1] < tspan[0]:
        sign = -1.0

    # Initialize h, the step size
    h = sign * abs(dt0)

    # Produce the initial X
    for i in range(0, statesize):
        X[i] = y0[i]

    v = funcptr(time, X)
    for i in range(0, statesize):
        dX[i] = v[i]

    for t in range(1, tspansize):
        goalt = tspan[t]
        while (sign * time) < (sign * goalt):
            if (sign * (time + h)) > (sign * goalt):
                h = goalt - time

            # Compute the slopes
            for i in range(0, statesize):
                f[i] = dX[i]
            for j in range(0, 6):
                for i in range(0, statesize):
                    Xcurrent[i] = 0.0
                    for k in range(0, 6):
                        Xcurrent[i] += f[i + k * statesize] * beta[k][j]
                    Xcurrent[i] = Xcurrent[i] * h + X[i]
                v = funcptr(time + alpha[j] * h, Xcurrent)
                for i in range(0, statesize):
                    f[i + (j + 1) * statesize] = v[i]

            for i in range(0, statesize):
                Xcurrent[i] = 0.0
                Xcurrent5[i] = 0.0
                for j in range(0, 6):
                    Xcurrent[i] += f[i + j * statesize] * chi4[j]
                    Xcurrent5[i] += f[i + j * statesize] * chi5[j]

                X4[i] = X[i] + h * Xcurrent[i]
                X5[i] = X[i] + h * Xcurrent5[i]

            # Estimate the error and the acceptable error
            delta = 0.0
            for i in range(0, statesize):
                gamma1 = abs(X5[i] - X4[i])
                if gamma1 > delta:
                    delta = gamma1

            tau = 1.0
            for i in range(0, statesize):
                if abs(X[i] > tau):
                    tau = abs(X[i])
            tau = tau * tol

            # Update the solution only if the error is acceptable
            if delta <= tau:
                time = time + h
                htol = EPS * abs(time)
                for i in range(0, statesize):
                    X[i] = X5[i]
                    dX[i] = f[i + 6 * statesize]

            # Update the step size
            if delta < 1.0e-16:
                delta = 1.0e-16
            h = 0.9 * h * (tau / delta) ** 0.2  # 1/5 = 0.2
            """ if the step size is really small (comparing the step size to
                the whole time span times the integration tolerance) or if it is
                NaN, terminate the integration so it doesn't go on endlessly
            """
            if abs(h) <= htol:
                print(
                    "\n{:s}: Step size below machine epsilon at current time (h={:f}).\n".format(
                        namestr, h
                    )
                )
            elif h != h:
                print(
                    "\n{:s}: Step size is NaN. Please debug.\n".format(namestr)
                )
                t = tspansize - 1
                time = goalt

        # Add this state to Xout
        for i in range(0, statesize):
            Y[t][i] = X[i]
        T[t] = time

    return T, Y
