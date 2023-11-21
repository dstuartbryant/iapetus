# Linear Dynamic Systems

## State-Space Representation
Denote the state of an object of interest at some time $t$ as a vector $\vec{x}(t)$. 

The unforced dynamics of a linear system can be represented by the following homogeneous matrix differential equation:
$$
\begin{equation}
\dot{\vec{x}}(t) = F(t)\vec{x}(t),
\end{equation}
$$

where, if $\vec{x}(t)$ is an $N\times 1$ vector, $F(t)$ is an $N\times N$ matrix which is sometimes referred to as the partial-derivative matrix [@vallado2013fundamentals-of].

## State Transition Matrix
Given a state $\vec{x}(t_0)$ at some initial time $t_0$ and assuming no forcing functions are applied, i.e., system inputs, the state $\vec{x}(t)$ at some other time $t$ can be computed via the *state transition matrix* $\Phi(t,t_0)$ via

$$
\begin{equation}
\vec{x}(t) = \Phi(t, t_0)\vec{x}_0.
\end{equation}
$$

The state transition matrix is defined by

$$
\begin{equation}
\dot{\Phi}(t, t_0) = F(t)\Phi(t,t_0), 
\end{equation}
$$

along with the following identities:

\begin{align}
\Phi(t_k,t_k) &= I, \\
\Phi(t_i,t_k) &= \Phi(t_i,t_j)\Phi(t_j,t_k), \\
\Phi(t_i,t_k) &= \Phi^{-1}(t_i,t_k).
\end{align}

### Solving for the State Transition Matrix
There are a variety of methods for solving for the state transition matrix, some of which are noted and/or discussed below.

#### Direct Differentiation
If the $F(t)$ matrix contains constant coefficients, then this approach is a good option.

The steps are 

1. Write the expressions for the system in state-space form to define $\dot{\vec{x}}(t)$, $\vec{x}(t)$, and $F(t)$, and define initial conditions,
2. Solve the system of equations through integration and applying intial conditions, resulting in fully formed expressions (equations) for each state vector component,
3. Then, the state transition matrix can be determined from direct differentiation as

\begin{equation}
\Phi(t,t_0) = \frac{\partial \vec{x}(t)}{\partial \vec{x}(t_0)}.
\end{equation}

#### Integrate STM Derivative Directly
With an $N\times 1$ state vector, the expression defining $\Phi$ can be expanded out

\begin{equation}
\begin{bmatrix}
\dot{\phi}_{11} & \dot{\phi}_{12} & \ldots &\dot{\phi}_{1N}\\
\dot{\phi}_{21} &  & \ldots & \vdots \\
\vdots & \ddots & & \\
\dot{\phi}_{N1} & \dot{\phi}_{N2} & \ldots &\dot{\phi}_{NN}\\
\end{bmatrix}
=
\begin{bmatrix}
F_{11} & F_{12} & \ldots &F_{1N}\\
F_{21} &  & \ldots & \vdots \\
\vdots & \ddots & & \\
F_{N1} & F_{N2} & \ldots &F_{NN}\\
\end{bmatrix}
\begin{bmatrix}
\phi_{11} & \phi_{12} & \ldots &\phi_{1N}\\
\phi_{21} &  & \ldots & \vdots \\
\vdots & \ddots & & \\
\phi_{N1} & \phi_{N2} & \ldots &\phi_{NN}\\
\end{bmatrix},
\end{equation}

and just carry out the matrix multiplication and solve for each $\phi$ by inegrating directly.

#### Laplace Tranforms Solution
I hate Laplace Transforms, so, see Born's Statisitical OD[@born] and Ogata's[@ogata] texts for more information on that.

#### Matrix Exponential for LTI Systems
Some systems are designed to be **time-invariant**, which means that when we
propagate a state through time, every time interval is identical. Linear systems such as these are sometimes referred to as linear time-invariant (LTI) systems.

For LTIs, the state transition matrix is 

\begin{equation}
\Phi(t,t_0) = \Phi(t-t_0),
\end{equation}

and its solution is somewhat simplified by leveraging matrix exponentials [@gelb] via

\begin{equation}
\Phi(t-t_0) = e^{F(t-t_0)},
\end{equation}

where the matrix exponetial of a matrix $A$ is defined as

\begin{equation}
e^A = I + A + \frac{A^2}{2!} + \frac{A^3}{3!} + \ldots
\end{equation}

#### Numerical Integration
...ADD NOTES HERE...


## Examples

### Straight Line Particle at Constant Velocity

Say that an object of interest is a particle moving along a straight line at constant velocity. The particle's position and velocity at some time $t$ are denoted as $r(t)$ and $v(t)$, respectively.

The state variables can be defined as

\begin{align}
x_1(t) &= r(t), \\
x_2(t) &= \dot{x}_1(t) = v(t).
\end{align}

Noting that $\dot{x}_2=0$ due to having zero acceleration, the dynamics of such a system are then represented as

\begin{equation}
\begin{bmatrix}
\dot{x}_1 \\
\dot{x}_2
\end{bmatrix}
=
\begin{bmatrix}
0 & 1\\
0 & 0
\end{bmatrix}
\begin{bmatrix}
x_1 \\
x_2
\end{bmatrix}.
\end{equation}

Alternatively, the $F(t)$ matrix can derived as $F(t) = \frac{\partial \dot{\vec{x}}(t)}{\partial \vec{x}(t)}$. So, for the straight line constant velocity example, 

\begin{equation}
F(t) = \frac{\partial \dot{\vec{x}}(t)}{\partial \vec{x}(t)} = 
\begin{bmatrix}
\frac{\partial \dot{x}_1}{\partial x_1} & \frac{\partial \dot{x}_1}{\partial x_2} \\
\frac{\partial \dot{x}_2}{\partial x_1} & \frac{\partial \dot{x}_2}{\partial x_2}
\end{bmatrix}
=
\begin{bmatrix}
0 & 1\\
0 & 0
\end{bmatrix}
\end{equation}

(Note that $\frac{\partial \dot{x}_1}{\partial x_1} =0$ due to the fact that we've defined the system such that velocity is constant, and therefore does not depend on, i.e., vary with respect with, the particle's position at any given time.)

Solving for $F(t)$ in this manner is why it is sometimes referred to as the ***Jacobian matrix***. 

Note, as per the definitions of this system, e.g., constant velocity, this system is an LTI, so solve for the state transition matrix accordingly.

\begin{equation}
\Phi(t-t_0) = e^{F(t-t_0)} = I + F(t-t_0) + \frac{(F(t-t_0))^2}{2!} + + \frac{(F(t-t_0))^3}{3!} + \ldots,
\end{equation}

but if we ignore the higher order terms, this reduces simply to,


\begin{align}
\Phi(t-t_0) &= 
\begin{bmatrix}
1 & 0 \\
0 & 1
\end{bmatrix}
+ 
\begin{bmatrix}
0 & 1 \\
0 & 0
\end{bmatrix}
\cdot (t-t_0) \\
&= 
\begin{bmatrix}
1 & (t-t_0) \\
0 & 1
\end{bmatrix}.
\end{align}

Sometimes the notation $\Delta t$ is used to denote a quantity such as $t-t_0$, and it's common to see the expression for $\Phi$ above as

\begin{equation}
\Phi = 
\begin{bmatrix}
1 & \Delta t \\
0 & 1
\end{bmatrix}.
\end{equation}

Given a state vector $\vec{x}(t_i) = \left[r_i, v_i \right]^T$ at time $t_i$, then the state vector $\vec{x}(t_j) = \left[r_j, v_j \right]^T$ at time $t_j$ can be computed, defining $\Delta t = t_j - t_i$, as

\begin{equation}
\begin{bmatrix}
r_j \\
v_j
\end{bmatrix}
 = 
\begin{bmatrix}
1 & \Delta t \\
0 & 1
\end{bmatrix}
\begin{bmatrix}
r_i \\
v_i
\end{bmatrix},
\end{equation}

\begin{align}
r_j &= r_i + v_i\Delta t, \\
v_j &= v_i.
\end{align}
