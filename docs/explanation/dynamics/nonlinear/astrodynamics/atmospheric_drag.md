# Atmospheric Drag

## Accelerations
Assuming a satellite is orbiting Earth, denote the satellite's state vector as

\begin{equation}
\vec{X} = [p_i, p_j, p_k, v_i, v_j, v_k]^T,
\end{equation}

where $p$ and $v$ denote position and velocity, respectively, and where $i$, $j$, and $k$ denote the axes of an Earth-centered inertial reference frame.

The acceleration due to atmospheric drag is defined as

$$
\label{eq:accel_formal} \tag{1}
\vec{a} = -\frac{1}{2}\frac{C_D A}{m}\rho v_{rel}\vec{v}_{rel}
$$

where $C_D$ is a satellite's drag coefficient, $A$ is the cross-sectional area of a satellite that is oriented in the direction of the satellite's velocity vector relative to the atmosphere, $m$ is the satellite's mass, $\rho$ is atmospheric density, and here we define $\vec{v}_{rel}$, the satellite's velocity relative to the atmosphere, as

\begin{equation} \label{eq:v_rel_with_derivatives} \tag{2}
\vec{v}_{rel} = \left[\frac{\mathrm{d}p_i}{\mathrm{d}t}+\omega_\oplus p_j, \frac{\mathrm{d}p_j}{\mathrm{d}t}-\omega_\oplus p_i, \frac{\mathrm{d}p_k}{\mathrm{d}t}\right]^T,
\end{equation}

where $\omega_\oplus$ denotes the Earth's rotational rate.

Sometimes it's useful to group all of the constant terms in Eq. \eqref{eq:accel_formal} into a single parameter called the ballistic coefficient, $B^*$, like so [@palmer2021]

\begin{equation} \label{eq:bstar} \tag{24}
B^* = \frac{1}{2}\frac{C_D A}{m},
\end{equation}

so then \eqref{eq:accel_formal} can be re-written as

$$
\vec{a} = -B^*\rho v_{rel}\vec{v}_{rel}. \label{eq:accel_formal_2} \tag{25}
$$

Since

$$
\frac{\mathrm{d}p_i}{\mathrm{d}t} = v_i, \;\; \frac{\mathrm{d}p_j}{\mathrm{d}t} = v_j,\;\;\frac{\mathrm{d}p_k}{\mathrm{d}t} =v_k,
$$

then \eqref{eq:v_rel_with_derivatives} can be simplified as

\begin{align}
\vec{v}_{rel} &= \left[v_i+\omega_\oplus p_j, v_j-\omega_\oplus p_i, v_k\right]^T,
\\
&= \left[ v_{rel,i}, v_{rel,j}, v_{rel,k}  \right] .\label{eq:rel_vel_vector} \tag{3}
\end{align}

Since 

\begin{equation} \label{eq:rel_vel_mag} \tag{4}
v_{rel} = |\vec{v}_{rel}| =\left( v_{rel,i}^2 + v_{rel,j}^2 + v_{rel,k}^2\right)^{1/2}, 
\end{equation}

\eqref{eq:accel_formal_2} can be re-written as

\begin{equation} \label{eq:accel_simplish} \tag{5}
\vec{a} = -B^*\rho \left( v_{rel,i}^2 + v_{rel,j}^2 + v_{rel,k}^2\right)^{1/2} \vec{v}_{rel}.
\end{equation}

Finally, note that atmospheric density as a function of a satellite's altitude, which makes it a function of it's position $p$, that is $\rho \equiv \rho(p)$, where

\begin{align}
\vec{p} &= \left[p_i, p_j, p_k \right]^T,
\\
p &= \left( p_i^2 + p_j^2 + p_k^2  \right)^{1/2},
\end{align}

and so, Eq. \eqref{eq:accel_simplish} can be re-written as

\begin{equation} \label{eq:accel_simplish_2} \tag{26}
\vec{a} = -B^*\rho(p) \left( v_{rel,i}^2 + v_{rel,j}^2 + v_{rel,k}^2\right)^{1/2} \vec{v}_{rel},
\end{equation}

which, given $\vec{a} = [a_i, a_j, a_k]^T$, can be broken out into components as

!!! note ""

    \begin{align}
    a_i &= -B^*\rho(p) \left( v_{rel,i}^2 + v_{rel,j}^2 + v_{rel,k}^2\right)^{1/2} v_{rel,i},
    \\
    a_j &= -B^*\rho(p) \left( v_{rel,i}^2 + v_{rel,j}^2 + v_{rel,k}^2\right)^{1/2} v_{rel,j},
    \\
    a_k &= -B^*\rho(p) \left( v_{rel,i}^2 + v_{rel,j}^2 + v_{rel,k}^2\right)^{1/2} v_{rel,k}.
    \end{align}

## Partial Derivatives

### General Atmospheric Density Model
For the purposes of these derivations, we only care about finding the partial derivatives of the acceleration components wrt all state variable components. See the appendix for the full derivations.

#### i-th component

!!! note ""

    \begin{align}
    \frac{\partial a_i}{\partial p_i} &= -B^* v_{rel,i} \left( \frac{\partial \rho}{\partial p_i} v_{rel} -   \frac{\omega_\oplus \rho v_{rel,j}}{v_{rel}}\right)
    \\
    \frac{\partial a_i}{\partial p_j} &= -B^*\left(\frac{\partial \rho}{\partial p_j} v_{rel} v_{rel,i} + \omega_\oplus\rho \left[ \frac{ v_{rel,i}^2}{v_{rel}} + v_{rel}\right]\right),
    \\
    \frac{\partial a_i}{\partial p_k} &= -B^*\frac{\partial \rho}{\partial p_k} v_{rel} v_{rel,i},
    \\
    \frac{\partial a_i}{\partial v_i} &= -B^*\rho\left( \frac{v_{rel,i}^2}{v_{rel}} + v_{rel}\right),
    \\
    \frac{\partial a_i}{\partial v_j} &= -B^*\rho \frac{v_{rel,i} v_{rel,j}}{v_{rel}},
    \\
    \frac{\partial a_i}{\partial v_k} &= -B^*\rho \frac{v_{rel,i} v_k}{v_{rel}},
    \end{align}

#### j-th component

!!! note ""

    \begin{align}
    \frac{\partial a_j}{\partial p_i} &= -B^*\left(\frac{\partial \rho}{\partial p_i} v_{rel} v_{rel,j} - \omega_\oplus\rho \left[ \frac{ v_{rel,j}^2}{v_{rel}} + v_{rel}\right]\right),
    \\
    \frac{\partial a_j}{\partial p_j} &=  -B^*v_{rel,j}\left(\frac{\partial \rho}{\partial p_j} v_{rel}  +  \frac{\omega_\oplus \rho v_{rel,i}}{v_{rel}} \right),
    \\
    \frac{\partial a_j}{\partial p_k} &=  -B^*\frac{\partial \rho}{\partial p_k} v_{rel} v_{rel,j},
    \\
    \frac{\partial a_j}{\partial v_i} &=  -B^* \rho \frac{v_{rel,i} v_{rel,j}}{v_{rel}},
    \\
    \frac{\partial a_j}{\partial v_j} &=  -B^* \rho \left[ \frac{v_{rel,j}^2}{v_{rel}} + v_{rel}\right],
    \\
    \frac{\partial a_j}{\partial v_k} &=  -B^* \rho  \frac{v_{rel,j} v_k}{v_{rel}},
    \end{align}

#### k-th component

!!! note ""

    \begin{align}
    \frac{\partial a_k}{\partial p_i} &= -B^*v_k\left(\frac{\partial \rho}{\partial p_i} v_{rel} -  \frac{\omega_\oplus \rho v_{rel,j}}{v_{rel}} \right),
    \\
    \frac{\partial a_k}{\partial p_j} &= -B^*v_k\left(\frac{\partial \rho}{\partial p_j} v_{rel} +  \frac{\omega_\oplus \rho v_{rel,i}}{v_{rel}}\right),
    \\
    \frac{\partial a_k}{\partial p_k} &= -B^*\frac{\partial \rho}{\partial p_k} v_{rel} v_k,
    \\
    \frac{\partial a_k}{\partial v_i} &=  -B^* \rho  \frac{v_{rel,i} v_k}{v_{rel}},
    \\
    \frac{\partial a_k}{\partial v_j} &=  -B^* \rho  \frac{v_{rel,j} v_k}{v_{rel}}  ,
    \\
    \frac{\partial a_k}{\partial v_k} &=  -B^* \rho \left[ \frac{v_k^2}{v_{rel}} + v_{rel}\right],
    \end{align}

#### $B^*$ Partial Derivatives
In instances which $B^*$ is included in the state vector, the following partial derivatives are also needed.

!!! note ""

    \begin{align}
    \frac{\partial a_i}{\partial B^*} &= \rho(p) \left( v_{rel,i}^2 + v_{rel,j}^2 + v_{rel,k}^2\right)^{1/2} v_{rel,i},
    \\
    \frac{\partial a_j}{\partial B^*} &= \rho(p) \left( v_{rel,i}^2 + v_{rel,j}^2 + v_{rel,k}^2\right)^{1/2} v_{rel,j},
    \\
    \frac{\partial a_k}{\partial B^*} &=\rho(p) \left( v_{rel,i}^2 + v_{rel,j}^2 + v_{rel,k}^2\right)^{1/2} v_{rel,k}.
    \end{align}

### Exponential Atmospheric Density Model
The equation for atmospheric density for the exponential model is given as [@vallado2013fundamentals-of]

$$
\rho = \rho_0 \exp\left[ - \frac{h_{ellp}-h_0}{H} \right],
$$

where reference density $\rho_0$ is used with the reference altitude $h_0$, the actual altitude (above the ellipsoid) $h_{ellp}$, and the scale height $H$.

If this model is used for Earth's atmosphere, then the expression for $h_{ellp}$ becomes

$$
h_{ellp} = p - R_\oplus,
$$

where $R_\oplus$ is the Earth's equatorial radius. Since $\rho$ is a function of altitude (i.e., $h_{ellp}$), we can rewrite the expression as

$$
\rho(p) = \rho_0 \exp\left[ - \frac{p - R_\oplus -h_0}{H} \right].
$$

Note that $\rho_0$, $h_0$, and $H$ are all technically functions of $p$ too, but they all come from a lookup table, therefore, for the purposes of differentiation we treat them as constants like we treat $R_\oplus$.

Hence, the partial derivative of $\rho(p)$ with respect to $p$ is

$$
\frac{\partial \rho(p)}{\partial p} = -\frac{\rho_0}{H} \exp\left[ - \frac{p - R_\oplus -h_0}{H} \right].
$$

Letting $\rho(p) \equiv \rho$ to simplify notation, then the expression from Appendix C become

\begin{align}
\frac{\partial \rho}{\partial p_i} &= -\frac{\rho_0 p_i}{Hp} \exp\left[ - \frac{p - R_\oplus -h_0}{H} \right],
\\
\frac{\partial \rho}{\partial p_j} &= -\frac{\rho_0 p_j}{Hp} \exp\left[ - \frac{p - R_\oplus -h_0}{H} \right],
\\
\frac{\partial \rho}{\partial p_k} &= -\frac{\rho_0 p_k}{Hp} \exp\left[ - \frac{p - R_\oplus -h_0}{H} \right].
\end{align}

Plugging these values into the component-wise partial derivatives in the General Atmospheric Density Model case yields the following.

#### i-th component

!!! note ""

    \begin{align}
    \frac{\partial a_i}{\partial p_i} &= -B^* v_{rel,i} \left( -\frac{\rho_0 p_i}{Hp} \exp\left[ - \frac{p - R_\oplus -h_0}{H} \right] v_{rel} -   \frac{\omega_\oplus \rho v_{rel,j}}{v_{rel}}\right)
    \\
    \frac{\partial a_i}{\partial p_j} &= -B^*\left(-\frac{\rho_0 p_j}{Hp} \exp\left[ - \frac{p - R_\oplus -h_0}{H} \right] v_{rel} v_{rel,i} + \omega_\oplus\rho \left[ \frac{ v_{rel,i}^2}{v_{rel}} + v_{rel}\right]\right),
    \\
    \frac{\partial a_i}{\partial p_k} &= B^*\frac{\rho_0 p_k}{Hp} \exp\left[ - \frac{p - R_\oplus -h_0}{H} \right] v_{rel} v_{rel,i},
    \\
    \frac{\partial a_i}{\partial v_i} &= -B^*\rho\left( \frac{v_{rel,i}^2}{v_{rel}} + v_{rel}\right),
    \\
    \frac{\partial a_i}{\partial v_j} &= -B^*\rho \frac{v_{rel,i} v_{rel,j}}{v_{rel}},
    \\
    \frac{\partial a_i}{\partial v_k} &= -B^*\rho \frac{v_{rel,i} v_k}{v_{rel}},
    \end{align}

#### j-th component

!!! note ""

    \begin{align}
    \frac{\partial a_j}{\partial p_i} &= -B^*\left(-\frac{\rho_0 p_i}{Hp} \exp\left[ - \frac{p - R_\oplus -h_0}{H} \right] v_{rel} v_{rel,j} - \omega_\oplus\rho \left[ \frac{ v_{rel,j}^2}{v_{rel}} + v_{rel}\right]\right),
    \\
    \frac{\partial a_j}{\partial p_j} &=  -B^*v_{rel,j}\left(-\frac{\rho_0 p_j}{Hp} \exp\left[ - \frac{p - R_\oplus -h_0}{H} \right] v_{rel}  +  \frac{\omega_\oplus \rho v_{rel,i}}{v_{rel}} \right),
    \\
    \frac{\partial a_j}{\partial p_k} &=  B^*\frac{\rho_0 p_k}{Hp} \exp\left[ - \frac{p - R_\oplus -h_0}{H} \right] v_{rel} v_{rel,j},
    \\
    \frac{\partial a_j}{\partial v_i} &=  -B^* \rho \frac{v_{rel,i} v_{rel,j}}{v_{rel}},
    \\
    \frac{\partial a_j}{\partial v_j} &=  -B^* \rho \left[ \frac{v_{rel,j}^2}{v_{rel}} + v_{rel}\right],
    \\
    \frac{\partial a_j}{\partial v_k} &=  -B^* \rho  \frac{v_{rel,j} v_k}{v_{rel}},
    \end{align}

#### k-th component

!!! note ""

    \begin{align}
    \frac{\partial a_k}{\partial p_i} &= -B^*v_k\left(-\frac{\rho_0 p_i}{Hp} \exp\left[ - \frac{p - R_\oplus -h_0}{H} \right] v_{rel} -  \frac{\omega_\oplus \rho v_{rel,j}}{v_{rel}} \right),
    \\
    \frac{\partial a_k}{\partial p_j} &= -B^*v_k\left(-\frac{\rho_0 p_j}{Hp} \exp\left[ - \frac{p - R_\oplus -h_0}{H} \right] v_{rel} +  \frac{\omega_\oplus \rho v_{rel,i}}{v_{rel}}\right),
    \\
    \frac{\partial a_k}{\partial p_k} &= B^*\frac{\rho_0 p_k}{Hp} \exp\left[ - \frac{p - R_\oplus -h_0}{H} \right] v_{rel} v_k,
    \\
    \frac{\partial a_k}{\partial v_i} &=  -B^* \rho  \frac{v_{rel,i} v_k}{v_{rel}},
    \\
    \frac{\partial a_k}{\partial v_j} &=  -B^* \rho  \frac{v_{rel,j} v_k}{v_{rel}}  ,
    \\
    \frac{\partial a_k}{\partial v_k} &=  -B^* \rho \left[ \frac{v_k^2}{v_{rel}} + v_{rel}\right],
    \end{align}



## Appendices


### Appendix A: Relative Velocity Vector and Magnitude Partial Derivatives

First let's take the partial derivatives of Eqs. \eqref{eq:rel_vel_vector} and \eqref{eq:rel_vel_mag} with respect to all of the state vector variables.

#### Relative Velocity Vector Partials

<div class="grid cards" markdown>

-   \begin{align}
    \frac{\partial v_{rel,i}}{\partial p_i} &= 0
    \\
    \frac{\partial v_{rel,i}}{\partial p_j} &= \omega_\oplus
    \\
    \frac{\partial v_{rel,i}}{\partial p_k} &= 0
    \\
    \frac{\partial v_{rel,i}}{\partial v_i} &= 1
    \\
    \frac{\partial v_{rel,i}}{\partial v_j} &= 0
    \\
    \frac{\partial v_{rel,i}}{\partial v_k} &= 0
    \end{align}

    

-   \begin{align}
    \frac{\partial v_{rel,j}}{\partial p_i} &= -\omega_\oplus
    \\
    \frac{\partial v_{rel,j}}{\partial p_j} &= 0
    \\
    \frac{\partial v_{rel,j}}{\partial p_k} &= 0
    \\
    \frac{\partial v_{rel,j}}{\partial v_i} &= 0
    \\
    \frac{\partial v_{rel,j}}{\partial v_j} &= 1
    \\
    \frac{\partial v_{rel,j}}{\partial v_k} &= 0
    \end{align}

-   \begin{align}
    \frac{\partial v_k}{\partial p_i} &= 0
    \\
    \frac{\partial v_k}{\partial p_j} &= 0
    \\
    \frac{\partial v_k}{\partial p_k} &= 0
    \\
    \frac{\partial v_k}{\partial v_i} &=0
    \\
    \frac{\partial v_k}{\partial v_j} &=0
    \\
    \frac{\partial v_k}{\partial v_k} &=1
    \\
    \end{align}

</div>



#### Relative Velocity Magnitude Partials

Recall the chain rule where the derivative of $f(g(x))$ is $f'(g(x))g'(x)$, and note that 

\begin{equation}
\frac{\mathrm{d}}{\mathrm{d}x} x^2=2x,
\end{equation}

\begin{equation}
\frac{\mathrm{d}}{\mathrm{d}x} x^{1/2}=\frac{1}{2}x^{-1/2}.
\end{equation}

First, note that

<div class="grid cards" markdown>

-   \begin{align}
    \frac{\partial }{\partial p_i} v_{rel,i}^2 & = 0
    \\
    \frac{\partial }{\partial p_j} v_{rel,i}^2 & = 2\omega_\oplus v_{rel,i}
    \\
    \frac{\partial }{\partial p_k} v_{rel,i}^2 & = 0
    \\
    \frac{\partial }{\partial v_i} v_{rel,i}^2 & = 2v_{rel,i}
    \\
    \frac{\partial }{\partial v_j} v_{rel,i}^2 & = 0
    \\
    \frac{\partial }{\partial v_k} v_{rel,i}^2 & =0
    \end{align}

    

-   \begin{align}
    \frac{\partial }{\partial p_i} v_{rel,j}^2 & = -2 \omega_\oplus v_{rel,j}
    \\
    \frac{\partial }{\partial p_j} v_{rel,j}^2 & = 0
    \\
    \frac{\partial }{\partial p_k} v_{rel,j}^2 & =  0
    \\
    \frac{\partial }{\partial v_i} v_{rel,j}^2 & = 0
    \\
    \frac{\partial }{\partial v_j} v_{rel,j}^2 & = 2 v_{rel,j}
    \\
    \frac{\partial }{\partial v_k} v_{rel,j}^2 & = 0
    \end{align}

-   \begin{align}
    \frac{\partial }{\partial p_i} v_k^2 & = 0
    \\
    \frac{\partial }{\partial p_j} v_k^2 & = 0
    \\
    \frac{\partial }{\partial p_k} v_k^2 & = 0
    \\
    \frac{\partial }{\partial v_i} v_k^2 & = 0
    \\
    \frac{\partial }{\partial v_j} v_k^2 & = 0
    \\
    \frac{\partial }{\partial v_k} v_k^2 & = 2v_k
    \end{align}

</div>

Then, 

\begin{align}
\frac{\partial v_{rel}}{\partial p_i} &= \frac{-\omega_\oplus v_{rel,j}}{v_{rel}}
\\
\frac{\partial v_{rel}}{\partial p_j} &= \frac{\omega_\oplus v_{rel,i}}{v_{rel}}
\\
\frac{\partial v_{rel}}{\partial p_k} &= 0
\\
\frac{\partial v_{rel}}{\partial v_i} &= \frac{v_{rel,i}}{v_{rel}}
\\
\frac{\partial v_{rel}}{\partial v_j} &= \frac{v_{rel,j}}{v_{rel}}
\\
\frac{\partial v_{rel}}{\partial v_k} &= \frac{v_k}{v_{rel}}
\\
\end{align}

## Appendix B: Acceleration Component Partial Derivatives

Recall that if

\begin{equation*}
f(x) = g(x)h(x),
\end{equation*}

that

\begin{equation*}
f'(x) = g'(x)h(x) + g(x)h'(x).
\end{equation*}

#### i-th component

<!-- #############################  dai/dpi  ############################## -->

\begin{align}
\frac{\partial a_i}{\partial p_i} &= -B^*\frac{\partial }{\partial p_i}\left[\rho v_{rel} v_{rel,i} \right],
\\
&= -B^*\left(\frac{\partial \rho}{\partial p_i} v_{rel} v_{rel,i} + \rho \frac{\partial }{\partial p_i}\left[ v_{rel} v_{rel,i}\right]\right),
\\
&= -B^*\left(\frac{\partial \rho}{\partial p_i} v_{rel} v_{rel,i} + \rho \left[ \frac{\partial v_{rel}}{\partial p_i} v_{rel,i} + v_{rel}\frac{\partial v_{rel,i}}{\partial p_i}\right]\right),
\\
&= -B^* \left( \frac{\partial \rho}{\partial p_i} v_{rel} v_{rel,i} + \rho \left[ \frac{-\omega_\oplus v_{rel,j}}{v_{rel}} v_{rel,i} + v_{rel} \cdot 0     \right] \right)
\\
&= -B^* v_{rel,i} \left( \frac{\partial \rho}{\partial p_i} v_{rel} -   \frac{\omega_\oplus \rho v_{rel,j}}{v_{rel}}\right)
\end{align}

<!-- #############################  dai/dpj  ############################## -->

\begin{align}
\frac{\partial a_i}{\partial p_j} &= -B^*\left(\frac{\partial \rho}{\partial p_j} v_{rel} v_{rel,i} + \rho \left[ \frac{\partial v_{rel}}{\partial p_j} v_{rel,i} + v_{rel}\frac{\partial v_{rel,i}}{\partial p_j}\right]\right),
\\
&= -B^*\left(\frac{\partial \rho}{\partial p_j} v_{rel} v_{rel,i} + \omega_\oplus\rho \left[ \frac{ v_{rel,i}^2}{v_{rel}} + v_{rel}\right]\right),
\end{align}


<!-- #############################  dai/dpk  ############################## -->

\begin{align}
\frac{\partial a_i}{\partial p_k} &= -B^*\left(\frac{\partial \rho}{\partial p_k} v_{rel} v_{rel,i} + \rho \left[ \frac{\partial v_{rel}}{\partial p_k} v_{rel,i} + v_{rel}\frac{\partial v_{rel,i}}{\partial p_k}\right]\right),
\\
&= -B^*\frac{\partial \rho}{\partial p_k} v_{rel} v_{rel,i},
\end{align}

<!-- #############################  dai/dvi  ############################## -->

\begin{align}
\frac{\partial a_i}{\partial v_i} &= -B^*\left(\frac{\partial \rho}{\partial v_i} v_{rel} v_{rel,i} + \rho \left[ \frac{\partial v_{rel}}{\partial v_i} v_{rel,i} + v_{rel}\frac{\partial v_{rel,i}}{\partial v_i}\right]\right),
\\&= -B^*\rho\left( \frac{v_{rel,i}^2}{v_{rel}} + v_{rel}\right),
\end{align}

<!-- #############################  dai/dvj  ############################## -->

\begin{align}
\frac{\partial a_i}{\partial v_j} &= -B^*\left(\frac{\partial \rho}{\partial v_j} v_{rel} v_{rel,i} + \rho \left[ \frac{\partial v_{rel}}{\partial v_j} v_{rel,i} + v_{rel}\frac{\partial v_{rel,i}}{\partial v_j}\right]\right),
\\
&= -B^*\rho \frac{v_{rel,i} v_{rel,j}}{v_{rel}},
\end{align}

<!-- #############################  dai/dvk  ############################## -->

\begin{align}
\frac{\partial a_i}{\partial v_k} &= -B^*\left(\frac{\partial \rho}{\partial v_k} v_{rel} v_{rel,i} + \rho \left[ \frac{\partial v_{rel}}{\partial v_k} v_{rel,i} + v_{rel}\frac{\partial v_{rel,i}}{\partial v_k}\right]\right),
\\
&= -B^*\rho \frac{v_{rel,i} v_k}{v_{rel}},
\end{align}



#### j-th component


<!-- #############################  daj/dpi  ############################## -->

\begin{align}
\frac{\partial a_j}{\partial p_i} &= -B^*\frac{\partial }{\partial p_i}\left[\rho v_{rel} v_{rel,j} \right],
\\
&= -B^*\left(\frac{\partial \rho}{\partial p_i} v_{rel} v_{rel,j} + \rho \frac{\partial }{\partial p_i}\left[ v_{rel} v_{rel,j}\right]\right),
\\
&= -B^*\left(\frac{\partial \rho}{\partial p_i} v_{rel} v_{rel,j} + \rho \left[ \frac{\partial v_{rel}}{\partial p_i} v_{rel,j} + v_{rel}\frac{\partial v_{rel,j}}{\partial p_i}\right]\right),
\\
&= -B^*\left(\frac{\partial \rho}{\partial p_i} v_{rel} v_{rel,j} + \rho \left[ \frac{-\omega_\oplus v_{rel,j}}{v_{rel}} v_{rel,j} - v_{rel}\omega_\oplus\right]\right),
\\
&= -B^*\left(\frac{\partial \rho}{\partial p_i} v_{rel} v_{rel,j} - \omega_\oplus\rho \left[ \frac{ v_{rel,j}^2}{v_{rel}} + v_{rel}\right]\right),
\end{align}

<!-- #############################  daj/dpj  ############################## -->

\begin{align}
\frac{\partial a_j}{\partial p_j} &=  -B^*\left(\frac{\partial \rho}{\partial p_j} v_{rel} v_{rel,j} + \rho \left[ \frac{\partial v_{rel}}{\partial p_j} v_{rel,j} + v_{rel}\frac{\partial v_{rel,j}}{\partial p_j}\right]\right),
\\
&=  -B^*v_{rel,j}\left(\frac{\partial \rho}{\partial p_j} v_{rel}  +  \frac{\omega_\oplus \rho v_{rel,i}}{v_{rel}} \right),
\end{align}

<!-- #############################  daj/dpk  ############################## -->

\begin{align}
\frac{\partial a_j}{\partial p_k} &=  -B^*\left(\frac{\partial \rho}{\partial p_k} v_{rel} v_{rel,j} + \rho \left[ \frac{\partial v_{rel}}{\partial p_k} v_{rel,j} + v_{rel}\frac{\partial v_{rel,j}}{\partial p_k}\right]\right),
\\
&=  -B^*\frac{\partial \rho}{\partial p_k} v_{rel} v_{rel,j},
\end{align}

<!-- #############################  daj/dvi  ############################## -->

\begin{align}
\frac{\partial a_j}{\partial v_i} &=  -B^* \rho \left[ \frac{\partial v_{rel}}{\partial v_i} v_{rel,j} + v_{rel}\frac{\partial v_{rel,j}}{\partial v_i}\right],
\\
&=  -B^* \rho \frac{v_{rel,i} v_{rel,j}}{v_{rel}},
\end{align}


<!-- #############################  daj/dvj  ############################## -->

\begin{align}
\frac{\partial a_j}{\partial v_j} &=  -B^* \rho \left[ \frac{\partial v_{rel}}{\partial v_j} v_{rel,j} + v_{rel}\frac{\partial v_{rel,j}}{\partial v_j}\right],
\\
&=  -B^* \rho \left[ \frac{v_{rel,j}^2}{v_{rel}} + v_{rel}\right],
\end{align}

<!-- #############################  daj/dvk  ############################## -->

\begin{align}
\frac{\partial a_j}{\partial v_k} &=  -B^* \rho \left[ \frac{\partial v_{rel}}{\partial v_k} v_{rel,j} + v_{rel}\frac{\partial v_{rel,j}}{\partial v_k}\right],
\\
&=  -B^* \rho  \frac{v_{rel,j} v_k}{v_{rel}},
\end{align}

#### k-th component

<!-- #############################  dak/dpi  ############################## -->

\begin{align}
\frac{\partial a_k}{\partial p_i} &= -B^*\frac{\partial }{\partial p_i}\left[\rho v_{rel} v_k \right],
\\
&= -B^*\left(\frac{\partial \rho}{\partial p_i} v_{rel} v_k + \rho \frac{\partial }{\partial p_i}\left[ v_{rel} v_k\right]\right),
\\
&= -B^*\left(\frac{\partial \rho}{\partial p_i} v_{rel} v_k + \rho \left[ \frac{\partial v_{rel}}{\partial p_i} v_k + v_{rel}\frac{\partial v_k}{\partial p_i}\right]\right),
\\
&= -B^*v_k\left(\frac{\partial \rho}{\partial p_i} v_{rel} -  \frac{\omega_\oplus \rho v_{rel,j}}{v_{rel}} \right),
\end{align}

<!-- #############################  dak/dpj  ############################## -->

\begin{align}
\frac{\partial a_k}{\partial p_j} &= -B^*\left(\frac{\partial \rho}{\partial p_j} v_{rel} v_k + \rho \left[ \frac{\partial v_{rel}}{\partial p_j} v_k + v_{rel}\frac{\partial v_k}{\partial p_j}\right]\right),
\\
&= -B^*v_k\left(\frac{\partial \rho}{\partial p_j} v_{rel} +  \frac{\omega_\oplus \rho v_{rel,i}}{v_{rel}}\right),
\end{align}

<!-- #############################  dak/dpk  ############################## -->

\begin{align}
\frac{\partial a_k}{\partial p_k} &= -B^*\left(\frac{\partial \rho}{\partial p_k} v_{rel} v_k + \rho \left[ \frac{\partial v_{rel}}{\partial p_k} v_k + v_{rel}\frac{\partial v_k}{\partial p_k}\right]\right),
\\
&= -B^*\frac{\partial \rho}{\partial p_k} v_{rel} v_k,
\end{align}

<!-- #############################  dak/dvi  ############################## -->

\begin{align}
\frac{\partial a_k}{\partial v_i} &=  -B^* \rho \left[ \frac{\partial v_{rel}}{\partial v_i} v_k + v_{rel}\frac{\partial v_k}{\partial v_i}\right],
\\
&=  -B^* \rho  \frac{v_{rel,i} v_k}{v_{rel}},
\end{align}

<!-- #############################  dak/dvj  ############################## -->

\begin{align}
\frac{\partial a_k}{\partial v_j} &=  -B^* \rho \left[ \frac{\partial v_{rel}}{\partial v_j} v_k + v_{rel}\frac{\partial v_k}{\partial v_j}\right],
\\
&=  -B^* \rho  \frac{v_{rel,j} v_k}{v_{rel}}  ,
\end{align}

<!-- #############################  dak/dvk  ############################## -->

\begin{align}
\frac{\partial a_k}{\partial v_k} &=  -B^* \rho \left[ \frac{\partial v_{rel}}{\partial v_k} v_k + v_{rel}\frac{\partial v_k}{\partial v_k}\right],
\\
&=  -B^* \rho \left[ \frac{v_k^2}{v_{rel}} + v_{rel}\right],
\end{align}

### Appendix C: Atmospheric Density Partial Derivatives

It can be handy to separate the partial derivative of atmospheric density like so [@palmer2021]

\begin{align}
\frac{\partial \rho}{\partial p_i} &= \frac{\partial \rho}{\partial p}\frac{\partial p}{\partial p_i} = \frac{\partial \rho}{\partial p}\frac{p_i}{p},
\\
\frac{\partial \rho}{\partial p_j} &= \frac{\partial \rho}{\partial p}\frac{\partial p}{\partial p_j} = \frac{\partial \rho}{\partial p}\frac{p_j}{p},
\\
\frac{\partial \rho}{\partial p_k} &= \frac{\partial \rho}{\partial p}\frac{\partial p}{\partial p_k} = \frac{\partial \rho}{\partial p}\frac{p_k}{p},
\end{align}

where $\rho(p) \equiv \rho$ for notational simplicity.
