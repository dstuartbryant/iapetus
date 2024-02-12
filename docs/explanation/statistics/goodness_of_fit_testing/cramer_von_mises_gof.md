# Cramer von Mises

## Background

Given a random sample $X_1, \ldots, X_n$ of size $n$, let (and/or ensure) $X_1 < X_2 < \ldots < X_n$, and further suppose that the distribution of $X$ is given by $F(x)$[@d2017goodness].

Denote the empirical distribution function (EDF) as $F_n(x)$ defined as[@d2017goodness]

$$
\begin{equation}
F_n(x) = \frac{\text{number of observations $\leq x$}}{n}; \;\; -\infty < x < \infty.
\end{equation}
$$

Let $\bar{F}(x; \theta)$ denote a continuous distribution where $\theta$ is a vector of parameters, e.g., for a normal distribution under test, $\theta = \left( \mu, \sigma \right)$. We're interested in testing if $\bar{F}(x; \theta)$ is a good approximation of our sample's true distribution, i.e., we want to know if $F(x) \approx \bar{F}(x;\theta)$. Furthermore, assume that $\theta$ is fully specified.


$$
\begin{equation} \label{eq:w2} \tag{1}
W^2 = \sum_i^n \left( Z_i - \frac{2i-1}{2n}   \right)^2 + \frac{1}{12n}
\end{equation}
$$

$$
\begin{equation} \label{eq:w2mod} \tag{2}
\hat{W}^2 = \left(W^2 - \frac{0.4}{n} + \frac{0.6}{n^2} \right) \left(1 + \frac{1}{n}  \right)
\end{equation}
$$

## Steps

These steps come from D’Agostino and Stephens [@d2017goodness], Section 4.4, pg. 104.

Given

* a random sample $X_1, \ldots, X_n$ of size $n$,
* a fully specified distribution to test against $\bar{F}(x; \theta)$,
* a chosen significance level $\alpha$,

Do 

1. Ensure the $X_i$ of our random sample are in ascending order, $X_1 < X_2 < \ldots < X_n$,
2. Calculate $Z_i = \bar{F}(X_i;\theta), i=1,\ldots,n$,
3. Compute $W^2$ test statistic via Eq. \eqref{eq:w2}, which is Eq. 4.2 from Reference 1 [@d2017goodness], using $Z$ values from Step 2 as input,
4. Compute the modified version of the test statistic using Eq. \eqref{eq:w2mod},
5. Consult Table 1 from Csörgő et al.[@csorgHo1996exact] 
    - find the appropriate row of the table based on $n$,
    - find the appropriate column of the table based on $\alpha = 1-p$, where $p$ denotes the p-value in the column header,
6. If the computed modified test statistic exceeds the value found in the row and column of the table, then the null hypothesis $H_0$ is rejected at significance level $\alpha$.



Ref 2 [@csorgHo1996exact]
Ref 3 [@stephens1970use]



