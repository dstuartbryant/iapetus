# Discrete Kalman Filter

## Nomenclature

* $n$, number of elements in state $\vec{X}$,
* $m$, number of elements in observation $\vec{z}$,
* $\hat{X}_k$, best estimate of state $\vec{X}$ at time $t_k$
* $\hat{P}_k$, best estimate of $n\times n$ covariance $P$ at time $t_k$
* $\hat{X}_{k+1|k}$, predicted state $\vec{X}$ at time $t_{k+1}$ given its estimate from time $t_k$,
* $\hat{P}_{k+1|k}$, predicted  $n\times n$ covariance $P$ at time $t_{k+1}$ given its estimate from time $t_k$,
* $Q_k$,  $n\times n$ process noise matrix at time $t_k$,
* $\Phi\left(t_{k+1}, t_k\right)$,  $n\times n$ state transition matrix (STM) from $t_k$ to $t_{k+1}$,
* $\vec{z}_{k+1}$, observation vector at time $t_{k+1}$,
* $R_{k+1}$, $m\times m$ observation covariance matrix at time $t_{k+1}$,
* $H_{k+1}$, $m\times n$ partial derivative matrix mapping state space to observation space at time $t_{k+1}$,
* $K_{k+1}$, $n \times m$ Kalman gain matrix at time $t_{k+1}$,
* $I_{n\times n}$, $n \times n$ identity matrix,
* $\hat{X}_{k+1}$, $n\times 1$ updated state $\vec{X}$ at time $t_{k+1}$,
* $P_{k+1}$, $n\times n$ updated covariance matrix $P$ at time $t_{k+1}$.

## Application

The *discrete* Kalman filter is applied to linear systems, for which some simplifications occur.

In nonlinear cases, $\Phi$ is better referred to as the *error* state transition matrix, as it is used to propagate errors of the estimate of state $\vec{X}$ through time, however, in the linear case, we regard $\Phi$ simply as the state transition matrix, as it can propagate the state exactly [@vallado2013fundamentals-of]. I.e., Taylor series expansion about the reference trajectory is unnecessary in this case and we use state $\vec{X}$ instead of a vector representing corrections to the state, e.g., $\delta \vec{x}$.

Specifics on process noise $Q$ are omitted here, but more detail on them can be found in ... `add reference here`

!!! bug

    Add reference to process noise in paragraph above.

There's also no use of a residual matrix $\tilde{b}$ as in the non-linear case. Instead, we simply use $\vec{z} - H\vec{X}$.

Additionally, the state is updated  each time a new observation is received and the $\Phi$ and $H$ matrices are almost always calculated  analytically.

For more details, see [@vallado2013fundamentals-of].

## Algorithm [@vallado2013fundamentals-of]

At time $t_{k+1}$ we have

* $\hat{X}_k$,
* $P_k$,
* $Q_k$,
* $\vec{z}_{k+1}$,
* $R_{k+1}$.

### Prediction

The state transition matrix $\Phi\left(t_{k+1}, t_k\right)$ is usually previously determined analytically via the relationship

$$
\Phi\left(t_{k+1}, t_k\right) = \frac{\partial \hat{X}_{k+1}}{\partial \hat{X}_k}
$$

and can typically be derived once and used for each filter iteration.

The predicted state and covariance are computed as

$$
\begin{align}
\hat{X}_{k+1|k} &= \Phi\left(t_{k+1}, t_k\right) \hat{X}_k,
\\
P_{k+1|k} &= \Phi\left(t_{k+1}, t_k\right) P_k \Phi\left(t_{k+1}, t_k\right)^T + Q_k.
\end{align}
$$

### Update

Similar to the state transition matrix, the $H_{k+1}$ matrix is usually previously determined analytically via the relationship

$$
H_{k+1} = \frac{\partial \vec{z}_{k+1}}{\partial \hat{X}_{k+1|k}}
$$

and can typically be derived once and used for each filter iteration.

The Kalman gain and updated state and covariance are computed as

$$
\begin{align}
K_{k+1} &= P_{k+1|k} H_{k+1}^T \left[ H_{k+1} P_{k+1|k} H_{k+1}^T + R_{k+1}\right]^{-1},
\\
\hat{X}_{k+1} &= \hat{X}_{k+1|k} + K_{k+1} \left[  \vec{z}_{k+1}-H_{k+1}\hat{X}_{k+1|k}  \right],
\\
P_{k+1} &= \left[ I_{n\times n} - K_{k+1}H_{k+1}\right] P_{k+1|k}.
\end{align}
$$
