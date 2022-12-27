


??? note "Notation"
    $\mathbb{X}$ and $x$ denote the state space and a single-target state, respectively, such that $x\in\mathbb{X}$. A multi-target state $X$ is a set of single-target states, i.e., $X \subset \mathbb{X}$. In a Bayesian framework the multi-target state is modeled as a random finite set (RFS) $\Xi$.

Let $f_{\Xi}(X)$ denote the multi-target probability distribution of an RFS $\Xi$. As detailed by Mahler[@mahler:statistical-multisource-multitarget-information-fusion::2007_a], the probability hypothesis density (PHD) is analogous to the notion of expected value, illustrated as
$$
\begin{equation}
D_{\Xi}(x) \triangleq \mathrm{E}\left[ \delta_{\Xi}(x)\right] = \int \delta_{X}(x) \cdot f_{\Xi}(X)\delta X,
\end{equation}
$$
with
$$
\delta_{X}(x) \triangleq \sum_{w\in X}\delta_{w}(x),
$$
where $\delta_{w}(x)$ is the Dirac measure at $w$.

**The PHD is not a probability density.** It is a density function on the single-target state space $\mathbb{X}$. Other names for the PHD include intensity, intensity density, and first-moment density[@mahler:statistical-multisource-multitarget-information-fusion::2007_a;@vo_and_ma:the-gaussian-mixture-probability-hypothesis-density-filter:itsp:2006_a].

For any region $S\subseteq X$, the expected number of targets of $\Xi$ in that region is the integral of $D_{\Xi}$ in that region
$$
\int_S D_{\Xi}(x)dx = \mathrm{E}\left[ |S\cap\Xi | \right].
$$



