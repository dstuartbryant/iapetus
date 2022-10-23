First, take a look at the Single Target Bayes Filter. This will provide background on the unerlying notions, mathematics, and notation to help frame discussion of the Probability Hypothesis Density (PHD) and its filter implementation.

## Single Target Bayes Filter


Say that a target *had* state \(\mathbf{x}'\) at time step $k$, and we're interested in its state \(\mathbf{x}\) at time step $k+1.$ At time $k+1$ we have a collection of observations $Z^{k+1}=Z^k \cup \mathbf{z_{k+1}}$, where $Z^k:\mathbf{z}_1, \ldots, \mathbf{z}_k$ is a time sequence of collected observations at time step $k$.

Define

$$
\begin{equation}
f_{k+1|k}(\mathbf{x})
    \triangleq 
        \int 
            f_{k+1|k}(\mathbf{x}|\mathbf{x}')
            \cdot 
            f_{k|k}(\mathbf{x}'|Z^k)d\mathbf{x}',
% 
\label{eq1} \tag{1}
\end{equation}
$$


$$\begin{equation}
f_{k+1|k+1}(\mathbf{x})
    \triangleq
        \frac{
                f_{k+1}(\mathbf{z}_{k+1}|\mathbf{x})
                \cdot
                f_{k+1|k}(\mathbf{x}|Z^k)
            }
            {
                \int f_{k+1}(\mathbf{z}_{k+1}|\mathbf{y})
                \cdot
                f_{k+1|k}(\mathbf{y}|Z^k)d\mathbf{y}
            },
% 
\label{eq2} \tag{2}
\end{equation}$$

where

- $f_{k|k}(\mathbf{x}'|Z^k)$ is the prior probablity density for the target was in state $\mathbf{x}'$,
- the probability density resulting from the product $f_{k+1|k}(\mathbf{x}|\mathbf{x}')\cdot f_{k|k}(\mathbf{x}'|Z^k)$ represents the target having state $\mathbf{x}$ at time $k+1$ given that it had state $\mathbf{x}'$ at time $k$,
- $\ldots$ to be continuted $\ldots$