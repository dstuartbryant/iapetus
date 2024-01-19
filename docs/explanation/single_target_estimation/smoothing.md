---
caption:
  table:
    position: "top"
---
# Smoothing

NEED TO ADD SMOOTHED COVARIANCE TOO.

Smoothing operation when using a sequential filter.

Let $\hat{x}_{k|\ell}$ denote the best estimate of $x$ at $t_k$ based on observations through $t_\ell$, where in general $\ell> k$.

Let

$$
\begin{align}
P_{k+1|k} &= \Phi(t_{k+1},t_k)P_{k|k} \Phi^T(t_{k+1},t_k) + \Gamma(t_{k+1},t_k)Q_k\Gamma^T(t_{k+1},t_k),
\\
S_k &= P_{k|k}\Phi^T(t_{k+1},t_k)\left( P_{k+1|k} \right)^{-1},
\end{align}
$$

then the smoothing algorithm is represented by

$$
\hat{x}_{k|\ell} = \hat{x}_{k|k} + S_k \left( \hat{x}_{k+1|\ell}  -  \Phi(t_{k+1},t_k)\hat{x}_{k|k} \right).
$$

Computation goes backward in index $k$, with $\hat{x}_{\ell|\ell}$, the filter
solution, as initial conditions.The following components of filter solutions when propagating "forward" must be stored so they can be used in smoothing computations:

* $\hat{x}_{k|k}$,
* $P_{k|k}$,
* $\Phi(t_{k+1},t_k)$,
* $\Gamma(t_{k+1},t_k)$,
    - note, that if $\Gamma$ is encapsulated in the version of $Q$ you're using for process noise, then you want to store all of the $Q_k$'s from each filter prediction step.

## Example with 3 Filter Iterations

!!! note

    In the tables below, $||$ used in notation because a single $|$ messes with the markdown table syntax.

Starting at time $t_0$ with initial conditions $\hat{x}_0$ and $P_0$, the elements in the right-most three columns of Table 1 below are stored.

Table: Forward Filter Iteration Components: Items in the 3 right-most columns are stored througout the filtering process to faciliate later smoothing.

|  $k$  |    $t_k$    |      $\hat{x}$        |       $P$     |  $\Phi$  |  $\Gamma$  |
| :---: | :-------: | :---------------------: |  :---------:  | :---------:  |  :---------:  |
|   1   |   $t_1$   |     $\hat{x}_{1\|1}$    |   $P_{1\|1}$  | $\Phi(t_{1},t_0)$ | $\Gamma(t_{1},t_0)$ |
|   2   |   $t_2$   |     $\hat{x}_{2\|2}$    |   $P_{2\|2}$  | $\Phi(t_{2},t_1)$ | $\Gamma(t_{2},t_1)$ |
|   3   |   $t_3$   |     $\hat{x}_{3\|3}$    |   $P_{3\|3}$  | $\Phi(t_{3},t_2)$ | $\Gamma(t_{3},t_2)$ |

After filtering "forward", we're at time $t_3$, i.e., $\ell=3$, and we set the following variables:

* $\hat{x}_{\ell|\ell} = \hat{x}_{3|3}$,
* $\hat{x}_{\ell-1|\ell-1} = \hat{x}_{2|2}$,
* $P_{\ell-1|\ell-1} = P_{2|2}$,
* $\Phi(t_\ell,t_{\ell-1}) = \Phi(t_3,t_2)$,
* $P_{\ell|\ell-1} = P_{3|2}$,

where

$$
P_{3|2} = \Phi(t_{3},t_2)P_{2|2} \Phi^T(t_{3},t_2) + \Gamma(t_{3},t_2)Q_2\Gamma^T(t_{3},t_2).
$$

Set $k=\ell-1$, i.e., $k=2$, and then compute

$$
\begin{align}
S_{\ell-1} = S_{2} &= P_{2|2}\Phi^T(t_{3},t_2)\left( P_{3|2} \right)^{-1},
\\
\hat{x}_{k|\ell} = \hat{x}_{\ell-1|\ell}= \hat{x}_{2|3} &= \hat{x}_{2|2} + S_2 \left( \hat{x}_{3|3}  -  \Phi(t_{3},t_2)\hat{x}_{2|2} \right).
\end{align}
$$

Now, from the "forward" filtering process we have the following variables:

* $\hat{x}_{\ell-2|\ell-2} = \hat{x}_{1|1}$,
* $P_{\ell-2|\ell-2} = P_{1|1}$,
* $\Phi(t_\ell-1,t_{\ell-2}) = \Phi(t_2,t_1)$,
* $P_{\ell-1|\ell-2} = P_{2|1}$,

where

$$
P_{2|1} = \Phi(t_{2},t_1)P_{1|1} \Phi^T(t_{2},t_1) + \Gamma(t_{2},t_1)Q_1\Gamma^T(t_{2},t_1),
$$

and, from the previous smoothing step, we have and use

* $\hat{x}_{\ell-1|\ell} = \hat{x}_{2|3}$.




Now set $k=\ell-2$, i.e., $k=1$, and compute

$$
\begin{align}
S_{\ell-2} = S_{1} &= P_{1|1}\Phi^T(t_2,t_1)\left( P_{2|1} \right)^{-1},
\\
\hat{x}_{k|\ell} = \hat{x}_{\ell-2|\ell}= \hat{x}_{1|3} &= \hat{x}_{1|1} + S_1 \left( \hat{x}_{2|3}  -  \Phi(t_{2},t_1)\hat{x}_{1|1} \right).
\end{align}
$$

And for this example, that's the last step. 

Now we have the best estimate of $x$ at time $t_1$ based on observations through times $t_1, t_2$ and $t_3$.

Notice that we don't iterate all the way back to $t_0$, because that wouldn't make sense due to the fact that $\hat{x}_0$ and $P_0$ are presumably not based on estimation, as they are our *a priori* information.

Therefore, we don't actually need to store (or use) $\Phi(t_{1},t_0)$ or $\Gamma(t_{1},t_0)$ - see Table 1. Nor do we use $P_{3|3}$ or $\Gamma(t_{3},t_2)$.

NEED TO ADD SMOOTHED COVARIANCE TOO.