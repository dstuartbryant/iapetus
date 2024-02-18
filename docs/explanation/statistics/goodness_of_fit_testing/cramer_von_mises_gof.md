---
caption:
  table:
    position: "top"
---
# Cramér von Mises

## Background

Given a random sample $X_1, \ldots, X_n$ of size $n$, let (and/or ensure) $X_1 < X_2 < \ldots < X_n$, and further suppose that the distribution of $X$ is given by $F(x)$[@d2017goodness].

Denote the empirical distribution function (EDF) as $F_n(x)$ defined as[@d2017goodness]

$$
\begin{equation}
F_n(x) = \frac{\text{number of observations $\leq x$}}{n}; \;\; -\infty < x < \infty.
\end{equation}
$$

Let $\bar{F}(x; \theta)$ denote a continuous distribution where $\theta$ is a vector of parameters, e.g., for a normal distribution under test, $\theta = \left( \mu, \sigma \right)$. We're interested in testing if $\bar{F}(x; \theta)$ is a good approximation of our sample's true distribution, i.e., we want to know if $F(x) \approx \bar{F}(x;\theta)$. 

There are different cases for which goodness-of-fit tests are performed. For the purposes of this document, we'll only focus on what is typically referred to as Case 0: the given distribution under test is fully specified [@d2017goodness], i.e., assume that the parameters encapsulated in $\theta$ are known/defined.

Expressions for the Cramér-von Mises statistic $W^2$ and the modified Cramér-von Mises statistic $\hat{W}^2$ are shown in Eqs. \eqref{eq:w2} and \eqref{eq:w2mod} below, respectively[@d2017goodness],

$$
\begin{equation} \label{eq:w2} \tag{1}
W^2 = \sum_i^n \left( Z_i - \frac{2i-1}{2n}   \right)^2 + \frac{1}{12n},
\end{equation}
$$

$$
\begin{equation} \label{eq:w2mod} \tag{2}
\hat{W}^2 = \left(W^2 - \frac{0.4}{n} + \frac{0.6}{n^2} \right) \left(1 + \frac{1}{n}  \right),
\end{equation}
$$

where the Z-values $Z_i$ are computed via the EDF of $X_i$, i.e., $F_n(x)$. Ensure that the Z-values are arranged in ascending order, i.e., $Z_1 < Z_2 < \ldots, < Z_n$ when using Eqs. \eqref{eq:w2} and \eqref{eq:w2mod}.

### Significance

Given a test statistic $T$ that takes the value $t$, the **significance level**, or **p-value**, of the statistic will then be the value $p = P(T >t)$, in the context of upper tail probability. The term can also be applied to the lower tail probability $P(T <t)$, however, in this context this is referred to as $q$, or **q-level**, i.e., $q=P(T <t)$, such that $q = 1-p$[@d2017goodness].

When a the value of a test statistic $t$ is less than its p-value $p = P(T >t)$ at a chosen significance level, we say that the test statistic is **not** significant (i.e., **insignificant**) and we **do not** reject the null hypothesis.

Presentation of these quantities varies from one reference to another, so caution must be taken to ensure appropriate statistical interpretations. 

For example, Table 4.2 from D’Agostino and Stephens [@d2017goodness] (pg. 105), contains a row of upper tail percentage points for the modified Cramér-von Mises statistic $\hat{W}^2$, arranged according to significance levels $\alpha$ - see Table 1 below. The significance levels are presented according to $p = P(T >t)$. In contrast, Table 1 from Csörgő et al.[@csorgHo1996exact] arranges percentage points for $\hat{W}^2$ according to, using their notation, values of $p$, such that entries in their Table 1 are $x$ such that $p = P\{w_n^2\leq x\}$. In the context of conventions applied in D’Agostino and Stephens regarding upper vs lower tail probabilities as discussed above, Csörgő's values of $p$ are actually values of $q$. Thus, if $\alpha = p$, and $q = 1-p$, then the columns of Csörgő's Table 1 are $1-\alpha$ - see Table 2 below.

!!! note "Assumed Rounding Error"

    The difference between the last column percentage point values in Tables 1 and 2 below, 1.167 and 1.6204, respectively, are assumed to be due to rounding *errors* between the two authors' numerical methods used when generating the tabulated data.

Table: Excerpt of Table 4.2 from D’Agostino and Stephens, Ref. [1], containing upper tail percentage points for modified Cramér-von Mises test statistic arranged by significance levels.

| Significance Level $\alpha$ | 0.25 | 0.15| 0.10 | 0.05| 0.025 | 0.01 | 0.001 |
| :---| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| Upper Taile Percentage Points | 0.209 | 0.284 | 0.347 | 0.461 | 0.581 | 0.743 | 1.167 |

!!! note "Table 2: $n$ > 1000"

    The percentage points in Table 2 are taken from the bottom row of Csörgő's Table 1, for sample size $n > 1000$.

Table: Excerpt of Table 1 from Csörgő et al., Ref. [2], containing percentage points for modified Cramér-von Mises test statistic arranged by q-values.

| Significance Level $1-\alpha$ | 0.75 | 0.85 | 0.90 | 0.95 | 0.975 | 0.99 | 0.999 |
| :---| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| Percentage Points | 0.20939  | 0.28406 | 0.34730 | 0.46136 | 0.58061 | 0.74346 | 1.16204 |


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
    - find the appropriate column of the table based on $\alpha = 1-p$, where $p$ denotes the p-value in the column header of Table 1,
6. If the computed modified test statistic exceeds the value found in the row and column of the table, then the null hypothesis $H_0$ is rejected at significance level $\alpha$.

## Example

```plotly
{"file_path": "explanation/statistics/goodness_of_fit_testing/figures/cramer_von_mises_recreate_fig_4_1.json"}
```


Ref 2 [@csorgHo1996exact]
Ref 3 [@stephens1970use]



