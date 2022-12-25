Let $f_{\Xi}(X)$ denote the multitarget probability distribution of a random finite set $\Xi$, where $X$ denotes some finite set (to include the empty set $\emptyset$) of single target state vectors $\mathbf{x}$, i.e., $\mathbf{x}\in X$. 

As detailed by Mahler[@mahler:statistical-multisource-multitarget-information-fusion::2007_a], the probability hypothesis density (PHD) is analogous to the notion of expected value, illustrated as
$$
\begin{equation}
D_{\Xi}(\mathbf{x}) \triangleq \mathrm{E}\left[ \delta_{\Xi}(\mathbf{x})\right] = \int \delta_{X}(\mathbf{x}) \cdot f_{\Xi}(X)\delta X,
\end{equation}
$$
with
$$
\delta_{X}(\mathbf{x}) \triangleq \sum_{\mathbf{w}\in X}\delta_{\mathbf{w}}(\mathbf{x}),
$$
where $\delta_{\mathbf{w}}(\mathbf{x})$ is the Dirac measure at $\mathbf{w}$.

??? note "Note on Dirac measure"
    To some, the presentation Mahler[@mahler:statistical-multisource-multitarget-information-fusion::2007_a] uses here regarding the Dirac measure *may* seem at odds with that used in other literature on the subject. In fact, he doesn't refer to $\delta_{\mathbf{w}}(\mathbf{x})$ as a Dirac measure, but as a *Dirac delta density*. Of course, presentations vary across the literature, but a salient example is taken from a reference Mahler often cites, Daley & Vere-Jones [@daley_and_vere-jone:an-introduction-to-the-theory-of-point-processes-volume-i:-elementary-theory-and-methods::2003_a], where they denote $\delta_\ell(\cdot)$ as the *Dirac measure* at $\ell$, being defined on Borel sets $A$ by 

    $$
    \delta_\ell(A) = 
    \begin{cases}
    1 & \text{if } x \in A, \\
    0 & \text{otherwise.}
    \end{cases}
    $$

    It may be possible to harmonize the two presentations if one considered the Borel $\sigma$-algebra $Y$ of the singleton $\{\mathbf{x}\}$.

See

- Daley, Vere-Jones
- REALLY, look at this! This [math stack exchange thread](https://math.stackexchange.com/questions/54197/can-a-dirac-delta-function-be-a-probability-density-function-of-a-random-variabl){:target="_blank"} seems appealing. REALLY, look at this!

[my-link]: http://google.com


<a href="http://google.com/" target="_blank">Hello, google!</a


