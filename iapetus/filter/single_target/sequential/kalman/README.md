# Kalman Filtering

## Linearized Kalman Filter (LKF)

*Notes from Stat OD 2 code.*

Looks like I start with initial state X0 at first observation epoch.

So at filter step k = 0 I have
* X_k = X_0
* P_k = P_0
* xhat_k = zero
* z_k
* R_k

and I don't perform a prediction step. I go straight into processing the first observation, but I don't save anything for state estimate information for the step k=0, except for the initial Q I used ,but likely b/c I forgot to add it in
outside of the function call to filter.

Next filter iteration step k = 1.

From the prediction step, I store for k = 1:

* X, the predicted X, as the state for X_k = X_1
* Phi, the STM,
* Gamma, the process noise STM,
* Q, process noise matrix,

then for k = 1 update step I store:

* xbar,
* P, 

after applying update equations.




k   X            Phi         Gamma        Q           xbar      P
-------------------------------------------------------------------------------
0   X0          I            N/A          Q0          0         P0
1   predicted   predicted    N/A          Q1 = Q0     updated   updated         
2   predicted   predicted    N/A          Q2 = Q0     updated   updated    

## Smoother

I flip all the output state info stashed from LKF.

Process the first step by itself, then iterate over remainig with for loop.