Say that a target *had* state \(\mathbf{x}'\) at time step $k$, and we're interested in its state \(\mathbf{x}\) at time step $k+1.$ At time $k+1$ we have a collection of observations $Z^{k+1}=Z^k \cup \mathbf{z_{k+1}}$, where $Z^k:\mathbf{z}_1, \ldots, \mathbf{z}_k$ is a time sequence of collected observations at time step $k$.

Define

$$
\begin{equation}
f_{k+1|k}(\mathbf{x}|Z^k)
    \triangleq 
        \int 
            f_{k+1|k}(\mathbf{x}|\mathbf{x}')
            \cdot 
            f_{k|k}(\mathbf{x}'|Z^k)d\mathbf{x}',
% 
\label{eq:sngl_tgt_bayes_predict} \tag{1}
\end{equation}
$$

$$\begin{equation}
f_{k+1}(\mathbf{z}_{k+1}|Z^k)
    \triangleq
        {
            \int f_{k+1}(\mathbf{z}_{k+1}|\mathbf{x})
            \cdot
            f_{k+1|k}(\mathbf{x}|Z^k)d\mathbf{x}
        },
% 
\label{eq:sngl_tgt_bayes_norm_factor} \tag{2}
\end{equation}$$



$$\begin{equation}
f_{k+1|k+1}(\mathbf{x}|Z^{k+1})
    \triangleq
        \frac{
                f_{k+1}(\mathbf{z}_{k+1}|\mathbf{x})
                \cdot
                f_{k+1|k}(\mathbf{x}|Z^k)
            }
            {
                f_{k+1}(\mathbf{z}_{k+1}|Z^k)
            },
% 
\label{eq:sngl_tgt_bayes_correct} \tag{3}
\end{equation}$$



where

- $f_{k|k}(\mathbf{x}'|Z^k)$ is the prior probablity density for *when* (i.e., $k$) the target was in state $\mathbf{x}'$ given the collection of observations $Z^k$ up to time $k$,
- the $f_{k+1|k}(\mathbf{x}|\mathbf{x}')$ pdf term encapsulates state transition (think *propagation*) process information, whereby the probability density resulting from the product $f_{k+1|k}(\mathbf{x}|\mathbf{x}')\cdot f_{k|k}(\mathbf{x}'|Z^k)$ represents the target having state $\mathbf{x}$ at time $k+1$ given that it had state $\mathbf{x}'$ at time $k$,


- the $f_{k+1}(\mathbf{z}_{k+1}|\mathbf{x})$ pdf term encapsulates the *likelihood* that observation $\mathbf{z}_{k+1}$ represents (i.e., *is associated with*) target $\mathbf{x}$ at time $k+1$.


Note that $\eqref{eq:sngl_tgt_bayes_correct}$ is constructed using Bayes' Rule where $\eqref{eq:sngl_tgt_bayes_norm_factor}$ is the Bayes normalization factor.