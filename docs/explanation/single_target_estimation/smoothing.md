---
caption:
  table:
    position: "top"
---
# Smoothing

Smoothing operation when using a sequential filter.

Let $\hat{x}_{k|\ell}$ denote the best estimate of $x$ at $t_k$ based on observations through $t_\ell$, where in general $\ell> k$.

Let

$$
\begin{align}
P_{k+1|k} &= \Phi(t_{k+1},t_k)P_{k|k} \Phi^T(t_{k+1},t_k) + \Gamma(t_{k+1},t_k)Q_k\Gamma^T(t_{k+1},t_k), \label{eq:int_cov} \tag{1}
\\
S_k &= P_{k|k}\Phi^T(t_{k+1},t_k)\left( P_{k+1|k} \right)^{-1}, \label{eq:int_innov} \tag{2}
\end{align}
$$

then the smoothing algorithm is represented by

$$
\begin{align}
\hat{x}_{k|\ell} &= \hat{x}_{k|k} + S_k \left( \hat{x}_{k+1|\ell}  -  \Phi(t_{k+1},t_k)\hat{x}_{k|k} \right), \label{eq:smoothed_x} \tag{3}
\\
P_{k|\ell} &= P_{k|k} + S_k \left( P_{k+1|\ell} - P_{k+1|k} \right)S_k^T. \label{eq:smoothed_cov} \tag{4}
\end{align}
$$

Computation goes backward in index $k$, with $\hat{x}_{\ell|\ell}$, the filter
solution, as initial conditions.The following components of filter solutions when propagating "forward" must be stored so they can be used in smoothing computations:

* $\hat{x}_{k|k}$,
* $P_{k|k}$,
* $\Phi(t_{k+1},t_k)$,
* $\Gamma(t_{k+1},t_k)$,
    - note, that if $\Gamma$ is encapsulated in the version of $Q$ you're using for process noise, then you want to store all of the $Q_k$'s from each filter prediction step.

## Example with 3 "Forward" Filter Iterations and 2 Smoothing Iterations

!!! note

    In the tables below, $||$ used in notation because a single $|$ messes with the markdown table syntax.

Starting at time $t_0$ with initial conditions $\hat{x}_0$ and $P_0$, the elements in the right-most three columns of Table 1 below are stored.

Table: Forward Filter Iteration Components: Items in the 3 right-most columns are stored througout the filtering process to faciliate later smoothing.

|  $k$  |    $t_k$    |      $\hat{x}$        |       $P$     |  $\Phi$  |  $\Gamma$  |
| :---: | :-------: | :---------------------: |  :---------:  | :---------:  |  :---------:  |
|   1   |   $t_1$   |     $\hat{x}_{1\|1}$    |   $P_{1\|1}$  | $\Phi(t_{1},t_0)$ | $\Gamma(t_{1},t_0)$ |
|   2   |   $t_2$   |     $\hat{x}_{2\|2}$    |   $P_{2\|2}$  | $\Phi(t_{2},t_1)$ | $\Gamma(t_{2},t_1)$ |
|   3   |   $t_3$   |     $\hat{x}_{3\|3}$    |   $P_{3\|3}$  | $\Phi(t_{3},t_2)$ | $\Gamma(t_{3},t_2)$ |

### First Smoothing Iteration

After filtering "forward", we're at time $t_3$, i.e., $\ell=3$, and when we begin the smoothing process, we set $k=\ell-1 = 2$ for the first iteration.

For Eq. \eqref{eq:int_cov} we have the following variables:

* $\Phi(t_{k+1},t_{k}) = \Phi(t_3,t_2)$,
* $P_{k|k} = P_{2|2}$,
* $\Gamma(t_{k+1},t_{k}) = \Gamma(t_3,t_2)$,
* and $Q_k=Q_2$.

And so using Eq. \eqref{eq:int_cov} we have

$$
P_{k+1|k} = P_{3|2} = \Phi(t_3,t_2)P_{2|2} \Phi^T(t_3,t_2) + \Gamma(t_3,t_2)Q_2\Gamma^T(t_3,t_2), \label{eq:P32} \tag{5}
$$

and then for Eq. \eqref{eq:int_innov} we have

$$
S_k = S_2 = P_{2|2} \Phi^T(t_3,t_2) \left( P_{3|2} \right)^{-1}. \label{eq:S2} \tag{6}
$$

For computing Eqs. \eqref{eq:smoothed_x} and \eqref{eq:smoothed_cov}, we have from "forward" filtering estimates

* $\hat{x}_{k+1|\ell}=\hat{x}_{3|3}$,
* $\hat{x}_{k|k}=\hat{x}_{2|2}$,
* $P_{k|k} = P_{2|2}$ (same as above for when computing Eqs. \eqref{eq:P32} and \eqref{eq:S2}),
* and $P_{k+1|\ell} = P_{3|3}$,

and using values from computing Eqs. \eqref{eq:P32} and \eqref{eq:S2}, we compute

$$
\begin{align}
\hat{x}_{k|\ell} = \hat{x}_{2|3}&= \hat{x}_{2|2} + S_2 \left( \hat{x}_{3|3}  -  \Phi(t_3,t_2)\hat{x}_{2|2} \right), \label{eq:smoothed_x1} \tag{7}
\\
P_{k|\ell} = P_{2|3} &= P_{2|2} + S_2 \left( P_{3|3} - P_{3|2} \right)S_2^T. \label{eq:smoothed_cov1} \tag{8}
\end{align}
$$

### Second Smoothing Iteration

Now we set $k = \ell-2 = 1$, and for computing Eq. \eqref{eq:int_cov} we have from the "forward" filtering estimates

* $\Phi(t_{k+1},t_{k}) = \Phi(t_2,t_1)$,
* $P_{k|k} = P_{1|1}$,
* $\Gamma(t_{k+1},t_{k}) = \Gamma(t_2,t_1)$,
* and $Q_k=Q_1$,

and then compute

$$
P_{k+1|k} = P_{2|1} = \Phi(t_2,t_1)P_{1|1} \Phi^T(t_2,t_1) + \Gamma(t_2,t_1)Q_1\Gamma^T(t_2,t_2). \label{eq:P21} \tag{9}
$$

For Eq. \eqref{eq:int_innov} we have

$$
S_k = S_1 = P_{1|1} \Phi^T(t_2,t_1) \left( P_{2|1} \right)^{-1}. \label{eq:S1} \tag{10}
$$

For computing Eqs. \eqref{eq:smoothed_x} and \eqref{eq:smoothed_cov}, we have from "forward" filtering estimates

* $\hat{x}_{k|k}=\hat{x}_{1|1}$,
* and $P_{k|k} = P_{1|1}$ (same as above for when computing Eqs. \eqref{eq:P21} and \eqref{eq:S1}).

Since this is the second iteration, we take the following values from the first smoothing iteration:

* $\hat{x}_{k+1|\ell}=\hat{x}_{2|3}$,
* and $P_{k+1|\ell} = P_{2|3}$.

Now, using the values from computing Eqs. \eqref{eq:P21} and \eqref{eq:S1}, we compute

$$
\begin{align}
\hat{x}_{k|\ell} = \hat{x}_{1|3}&= \hat{x}_{1|1} + S_1 \left( \hat{x}_{2|3}  -  \Phi(t_2,t_1)\hat{x}_{1|1} \right), \label{eq:smoothed_x2} \tag{11}
\\
P_{k|\ell} = P_{1|2} &= P_{1|1} + S_1 \left( P_{2|3} - P_{2|1} \right)S_1^T. \label{eq:smoothed_cov2} \tag{12}
\end{align}
$$

And for this example, that's the last step. 

Now we have the best estimate of $x$ at time $t_1$ based on observations through times $t_1, t_2$ and $t_3$.

Notice that we don't iterate all the way back to $t_0$, because that wouldn't make sense due to the fact that $\hat{x}_0$ and $P_0$ are presumably not based on estimation, as they are our *a priori* information.

Therefore, we don't actually need to store (or use) $\Phi(t_{1},t_0)$ or $\Gamma(t_{1},t_0)$ - see Table 1, however, it may make sense to store them and pass them into a smoothing algorithm if it makes index processing easier. 

