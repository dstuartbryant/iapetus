"""First attempt at utilizing sympy to generate LaTeX and Python content.

This is more of a scratch/draft file.

* Math to be done: Derive the partial derivaties matrix (10-41) on pg. 805 of
  Valldo (3rd Ed.).

Want to:
* generate usable LaTeX code for documentation,
* generate usable python code for mathmatical operations

"""

import sympy
from sympy import Matrix, cos, sin
from sympy.concrete.summations import summation

# from sympy.core.numbers import Infinity
from sympy.vector import CoordSys3D

# r_I, r_J, r_K, v_I, v_J, v_K, r_vec, v_vec, mu = sympy.symbols("r_I r_J r_K v_I v_J v_K r_vec v_vec mu")

x, y, z, v, u, w, mu = sympy.symbols("x y z v u w mu")

r = (x**2 + y**2 + z**2) ** (1 / 2)

ax = -mu * x / r**3
ay = -mu * y / r**3
az = -mu * z / r**3

da_dx = sympy.diff(ax, x)

# Convert to string to facilitate easier subsitituion
da_dx_str = str(da_dx)


# Make simplifying substitutions
def sub_r(string):
    r3 = ("(x**2 + y**2 + z**2)**1.5", "r**3")
    r5 = ("(x**2 + y**2 + z**2)**2.5", "r**5")
    if r3[0] in string:
        string = string.replace(r3[0], r3[1])
    if r5[0] in string:
        string = string.replace(r5[0], r5[1])
    return string


def substitute_r_magnitude(string):
    r_relations = [
        ("(x**2 + y**2 + z**2)**1.5", "r**3"),
        ("(x**2 + y**2 + z**2)**2.5", "r**5"),
    ]
    for rel in r_relations:
        if rel[0] in string:
            string = string.replace(rel[0], rel[1])
    return string


da_dx_str = sub_r(da_dx_str)

# Convert back to sympy object
da_dx_new = sympy.sympify(da_dx_str)

# Convert to LaTeX string
da_dx_latex = sympy.latex(da_dx_new)


# Attempt at taking Jacobian of matrix wrt. vector
X_dot = Matrix([u, v, w, ax, ay, az])
X = Matrix([x, y, z, u, v, w])
dXdot_dX = X_dot.jacobian(X)
dXdot_dX_str = str(dXdot_dX)
dXdot_dX_str = substitute_r_magnitude(dXdot_dX_str)
dXdot_dX = sympy.sympify(dXdot_dX_str)
dXdot_dX_latex = sympy.latex(dXdot_dX)


# Attempt at taking derivative of dU/dr with respect to r
(
    r,
    mu,
    R_Earth,
    P,
    C,
    S,
    phi,
    lamb,
    j,
    m,
    infinity,
) = sympy.symbols("r mu  R_Earth P C S phi lamb j m infinity")


inner_loop_fcn = (
    (R_Earth / r) ** j
    * (j + 1)
    * P
    * sin(phi)
    * (C * cos(m * lamb) + S * sin(m * lamb))
)

inner_loop_sum = summation(inner_loop_fcn, (m, 0, j))
outer_loop_sum = summation(inner_loop_sum, (j, 2, infinity))

dUdr = -mu / r**2 * outer_loop_sum

d2Udr2 = sympy.diff(dUdr, r)
