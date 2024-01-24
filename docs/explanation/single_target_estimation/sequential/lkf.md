# Linearized Kalman Filter

## Nomenclature

* $n$, number of elements in state $\vec{X}$,
* $m$, number of elements in observation $\vec{z}$,
* $\Delta t$, difference between to times, e.g., $\Delta t = t_{k+1} - t_k$,
* $\hat{X}_0$, $n\times 1$ *a priori* initial state vector at time $t_0$,
* $\delta\vec{x}_0$, $n\times 1$ initial state error vector at time $t_0$,
* $P_0$, $n\times n$ *a priori* initial covariance matrix at time $t_0$,
* $\delta\vec{x}_{k}$, $n\times 1$ state error vector at time $t_k$ where $k > 0$,
* $\hat{P}_k$, $n\times n$ covariance $P$ at time $t_k$,
* $\delta\vec{x}_{k+1|k}$, $n\times 1$ predicted state error vector at time $t_{k+1}$,
* $\hat{X}_{k+1|k}$, predicted state $\vec{X}$ at time $t_{k+1}$,
* $\hat{P}_{k+1|k}$, predicted  $n\times n$ covariance $P$ at time $t_{k+1}$ given its estimate from time $t_k$,
* $Q_k$,  $n\times n$ process noise matrix at time $t_k$,
* $\Phi\left(t_{k+1}, t_k\right)$,  $n\times n$ error state transition matrix (STM) from $t_k$ to $t_{k+1}$,
* $F$, $n\times n$ partial derivatives (Jacobian) matrix of "state rates" [@vallado2013fundamentals-of],
* $I_{n\times n}$, $n \times n$ identity matrix,
* $\vec{z}_{k+1}$, observation vector at time $t_{k+1}$,
* $R_{k+1}$, $m\times m$ observation covariance matrix at time $t_{k+1}$,
* $H_{k+1}$, $m\times n$ partial derivative matrix mapping state space to observation space at time $t_{k+1}$,
* $\tilde{b}_{k+1}$, $m\times 1$ residuals matrix,
* $K_{k+1}$, $n \times m$ Kalman gain matrix at time $t_{k+1}$,
* $\delta\vec{x}_{k+1}$, $n\times 1$ updated state error vector at time $t_{k+1}$,
* $\hat{X}_{k+1}$, $n\times 1$ updated state $\vec{X}$ at time $t_{k+1}$,
* $P_{k+1}$, $n\times n$ updated covariance matrix $P$ at time $t_{k+1}$.

## Application

The linear (i.e., *linearized*) Kalman filter (LKF) is an adaptation of the (discrete) Kalman filter  for nonlinear systems.

A key assumption is that the *a priori* state generated reference trajectory  is sufficiently close to the actual (true) trajectory.

The *a priori* state at time $t_0$, i.e., $\hat{X}_0$ and $P_0$, are used to define the reference trajectory across all timestamps for which observations are received. I.e., for all $t_k$ where $k>0$. This means that the reference trajectory is not recalculated at each filter iteration.

Since the initial state is used for all initial calculations, the state and the update steps are separated, but is also why the LKF cannot handle large changes in state like the Extended Kalman filter (EKF) can ... `insert reference here`

!!! bug

    Add link to EKF page.

!!! warning
    
    We do not use the updated state in subsequent predictions. All computations are based on a nominal trajectory without incorporating updates.

For the first iteration, the initial state vector error is assumed to be zero, i.e., $\delta\vec{x}_0 = 0_{n\times 1}$.

The error state transition matrix $\Phi$ can be found by numerical integration, however, for most LKF filters the following approximation of the Jacobian matrix can be used:

$$
\Phi\left(t_{k+1}, t_k\right) = I_{n\times n} + F(t_{k+1}-t_k) + \frac{1}{2!}F^2(t_{k+1}-t_k)^2 + \frac{1}{3!}F^3(t_{k+1}-t_k)^3 + \ldots.
$$

Specifics on process noise $Q$ are omitted here, but more detail on them can be found in ... `add reference here`

!!! bug

    Add reference to process noise in paragraph above.

The $\tilde{b}$ vector contains the residuals, and note that we use $H\delta{x}$ in this implementation, just as in nonlinear least squares ... `add reference here`

!!! bug

    Add reference to nonlinear least squares in paragraph above.

### Good For

* systems that require current estimates of the state and computational speed,
* autonomous navigation systems and similar tasks,

### Draw Backs

* errors will grow with time (the reference trajectory will diverge from the truth) and at some point the state must be reinitialized and a new reference trajectory must be generated

## Algorithm [@vallado2013fundamentals-of]

At time $t_0$ we have $\delta\vec{x}_0$, $\hat{X}_0$, $P_0$ and $Q_0$.

Determine the Jacobian matrix $F$ once at the outset. See discussion of the Jacobian matrix in [Linear Dynamic Systems](../../dynamics/linear/top_linear_dynamic_systems.md){:target="_blank"}.

At this stage we set $k=0$ such that 

$$
\begin{align}
\delta\vec{x}_k &= \delta\vec{x}_0 = 0_{n\times 1}, 
\\
\hat{X}_k &= \hat{X}_0,
\\
P_k &= P_0,
\\
Q_k &= Q_0,
\end{align}
$$

we have a collection of obserations from $t_{k+1}, \ldots, t_L$, wereh $L$ denotes the number of observations at different points in time, i.e., we have $\vec{z}_1, \ldots, \vec{z}_L$.

We can compute all of the predicted states and error state transition matrices from $t_{k+1}=t_1$ to $t_L$ given some chosen form of *propagator*. E.g., if we choose to integrate to find all $\Phi$s, we can form sets of pre-computed states and error state transition matrices


$$
\begin{align}
\mathbf{X} &= \left\{\hat{X}_{k+1|k} \; : \; k+1 = 1, \ldots, L  \right\},
\\
\mathbf{\Phi} &= \left\{\Phi\left(t_{k+1}, t_k\right) \; : \; k+1 = 1, \ldots, L  \right\}.
\end{align}
$$





* $Q_k$,
* $\vec{z}_{k+1}$,
* $R_{k+1}$.