import math

import numpy as np

EPS = np.spacing(1)


def rk78(funcptr, y0, tspan, dt, tspantol, extras=None):
    dt0 = dt
    tol = tspantol

    namestr = "rk78"

    statesize = len(y0)
    tspansize = len(tspan)
    Y = np.zeros((tspansize, statesize))
    T = np.zeros((tspansize,))
    f = np.zeros((statesize * 13,))
    X = np.zeros((statesize,))
    X7 = np.zeros((statesize,))
    X8 = np.zeros((statesize,))
    Xcurrent = np.zeros((statesize,))
    Xcurrent8 = np.zeros((statesize,))

    # Set coefficient matrices
    alpha = np.zeros((12,))
    alpha[0] = 1.0 / 18.0
    alpha[1] = 1.0 / 12.0
    alpha[2] = 1.0 / 8.0
    alpha[3] = 5.0 / 16.0
    alpha[4] = 3.0 / 8.0
    alpha[5] = 59.0 / 400.0
    alpha[6] = 93.0 / 200.0
    alpha[7] = 5490023248.0 / 9719169821.0
    alpha[8] = 13.0 / 20.0
    alpha[9] = 1201146811.0 / 1299019798.0
    alpha[10] = 1.0
    alpha[11] = 1.0

    beta = np.zeros((12, 12))
    beta[0][0] = 1.0 / 18.0
    beta[0][1] = 1.0 / 48.0
    beta[0][2] = 1.0 / 32.0
    beta[0][3] = 5.0 / 16.0
    beta[0][4] = 3.0 / 80.0
    beta[0][5] = 29443841.0 / 614563906.0
    beta[0][6] = 16016141.0 / 946692911.0
    beta[0][7] = 39632708.0 / 573591083.0
    beta[0][8] = 246121993.0 / 1340847787.0
    beta[0][9] = -1028468189.0 / 846180014.0
    beta[0][10] = 185892177.0 / 718116043.0
    beta[0][11] = 403863854.0 / 491063109.0

    beta[1][1] = 1.0 / 16.0

    beta[2][2] = 3.0 / 32.0
    beta[2][3] = -75.0 / 64.0

    beta[3][3] = 75.0 / 64.0
    beta[3][4] = 3.0 / 16.0
    beta[3][5] = 77736538.0 / 692538347.0
    beta[3][6] = 61564180.0 / 158732637.0
    beta[3][7] = -433636366.0 / 683701615.0
    beta[3][8] = -37695042795.0 / 15268766246.0
    beta[3][9] = 8478235783.0 / 508512852.0
    beta[3][10] = -3185094517.0 / 667107341.0
    beta[3][11] = -5068492393.0 / 434740067.0

    beta[4][4] = 3.0 / 20.0
    beta[4][5] = -28693883.0 / 1125000000.0
    beta[4][6] = 22789713.0 / 633445777.0
    beta[4][7] = -421739975.0 / 2616292301.0
    beta[4][8] = -309121744.0 / 1061227803.0
    beta[4][9] = 1311729495.0 / 1432422823.0
    beta[4][10] = -477755414.0 / 1098053517.0
    beta[4][11] = -411421997.0 / 543043805.0

    beta[5][5] = 23124283.0 / 1800000000.0
    beta[5][6] = 545815736.0 / 2771057229.0
    beta[5][7] = 100302831.0 / 723423059.0
    beta[5][8] = -12992083.0 / 490766935.0
    beta[5][9] = -10304129995.0 / 1701304382.0
    beta[5][10] = -703635378.0 / 230739211.0
    beta[5][11] = 652783627.0 / 914296604.0

    beta[6][6] = -180193667.0 / 1043307555.0
    beta[6][7] = 790204164.0 / 839813087.0
    beta[6][8] = 6005943493.0 / 2108947869.0
    beta[6][9] = -48777925059.0 / 3047939560.0
    beta[6][10] = 5731566787.0 / 1027545527.0
    beta[6][11] = 11173962825.0 / 925320556.0

    beta[7][7] = 800635310.0 / 3783071287.0
    beta[7][8] = 393006217.0 / 1396673457.0
    beta[7][9] = 15336726248.0 / 1032824649.0
    beta[7][10] = 5232866602.0 / 850066563.0
    beta[7][11] = -13158990841.0 / 6184727034.0

    beta[8][8] = 123872331.0 / 1001029789.0
    beta[8][9] = -45442868181.0 / 3398467696.0
    beta[8][10] = -4093664535.0 / 808688257.0
    beta[8][11] = 3936647629.0 / 1978049680.0

    beta[9][9] = 3065993473.0 / 597172653.0
    beta[9][10] = 3962137247.0 / 1805957418.0
    beta[9][11] = -160528059.0 / 685178525.0

    beta[10][10] = 65686358.0 / 487910083.0
    beta[10][11] = 248638103.0 / 1413531060.0

    chi7 = np.zeros((13,))
    chi7[0] = 13451932.0 / 455176623.0
    chi7[1] = 0.0
    chi7[2] = 0.0
    chi7[3] = 0.0
    chi7[4] = 0.0
    chi7[5] = -808719846.0 / 976000145.0
    chi7[6] = 1757004468.0 / 5645159321.0
    chi7[7] = 656045339.0 / 265891186.0
    chi7[8] = -3867574721.0 / 1518517206.0
    chi7[9] = 465885868.0 / 322736535.0
    chi7[10] = 53011238.0 / 667516719.0
    chi7[11] = 2.0 / 45.0
    chi7[12] = 0.0

    chi8 = np.zeros((13,))
    chi8[0] = 14005451.0 / 335480064.0
    chi8[1] = 0.0
    chi8[2] = 0.0
    chi8[3] = 0.0
    chi8[4] = 0.0
    chi8[5] = -59238493.0 / 1068277825.0
    chi8[6] = 181606767.0 / 758867731.0
    chi8[7] = 561292985.0 / 797845732.0
    chi8[8] = -1041891430.0 / 1371343529.0
    chi8[9] = 760417239.0 / 1151165299.0
    chi8[10] = 118820643.0 / 751138087.0
    chi8[11] = -528747749.0 / 2220607170.0
    chi8[12] = 1.0 / 4.0

    # Set initial points in the output variables
    for i in range(0, statesize):
        Y[0][i] = y0[i]
    T[0] = tspan[0]

    # Initialize f
    for i in range(0, statesize):
        for j in range(0, 13):
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

    for t in range(1, tspansize):
        goalt = tspan[t]
        while (sign * time) < (sign * goalt):
            if (sign * (time + h)) > (sign * goalt):
                h = goalt - time

            # Compute the slopes
            v = funcptr(time, X)
            for i in range(0, statesize):
                f[i] = v[i]
            for j in range(0, 12):
                for i in range(0, statesize):
                    Xcurrent[i] = 0.0
                    for k in range(0, 12):
                        Xcurrent[i] += f[i + k * statesize] * beta[k][j]
                    Xcurrent[i] = Xcurrent[i] * h + X[i]
                v = funcptr(time + alpha[j] * h, Xcurrent)
                for i in range(0, statesize):
                    f[i + (j + 1) * statesize] = v[i]

            for i in range(0, statesize):
                Xcurrent[i] = 0.0
                Xcurrent8[i] = 0.0
                for j in range(0, 13):
                    Xcurrent[i] += f[i + j * statesize] * chi7[j]
                    Xcurrent8[i] += f[i + j * statesize] * chi8[j]

                X7[i] = X[i] + h * Xcurrent[i]
                X8[i] = X[i] + h * Xcurrent8[i]

            # Estimate the error and the acceptable error
            delta = 0.0
            for i in range(0, statesize):
                gamma1 = abs(X8[i] - X7[i])
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
                    X[i] = X8[i]

            # Update the step size
            if delta < 1.0e-16:
                delta = 1.0e-16
            h = 0.9 * h * (tau / delta) ** 0.125  # 1/8 = 0.125
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
