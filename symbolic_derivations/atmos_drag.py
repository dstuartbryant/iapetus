"""Symbolic derivation of partial derivatives for satellite atmospheric drag
dynamics.
"""
import sympy

(
    p_i,
    p_j,
    p_k,
    v_i,
    v_j,
    v_k,
    R_E,
    w,
    Cd,
    A,
    m,
    rho_0,
    R_E,
    h0,
    H,
    r,
) = sympy.symbols(
    "p_i  p_j  p_k  v_i  v_j  v_k  R_E  w  Cd  A  m rho_0 R_E h0 H r"
)

# r = (p_i**2 + p_j**2 + p_k**2) ** (1 / 2)
rho = rho_0 * sympy.exp(-(r - R_E - h0) / H)

f = -1 / 2 * Cd * A * rho / m

v_rel = [v_i + w * p_j, v_j - w * p_i, v_k]
v_rel_mag = (v_rel[0] ** 2 + v_rel[1] ** 2 + v_rel[2] ** 2) ** (1 / 2)

ai = f * v_rel_mag**2 * v_rel[0] / v_rel_mag
aj = f * v_rel_mag**2 * v_rel[1] / v_rel_mag
ak = f * v_rel_mag**2 * v_rel[2] / v_rel_mag

drho_dp = sympy.diff(rho, r)
# drho_dpi = sympy.diff(rho, p_i),


# dai_dpi = sympy.diff(ai, p_i)
# dai_dpj = sympy.diff(ai, p_j)
# dai_dpk = sympy.diff(ai, p_k)
# dai_dvi = sympy.diff(ai, v_i)
# dai_dvj = sympy.diff(ai, v_j)
# dai_dvk = sympy.diff(ai, v_k)


# daj_dpi = sympy.diff(aj, p_i)
# daj_dpj = sympy.diff(aj, p_j)
# daj_dpk = sympy.diff(aj, p_k)
# daj_dvi = sympy.diff(aj, v_i)
# daj_dvj = sympy.diff(aj, v_j)
# daj_dvk = sympy.diff(aj, v_k)

# dak_dpi = sympy.diff(ak, p_i)
# dak_dpj = sympy.diff(ak, p_j)
# dak_dpk = sympy.diff(ak, p_k)
# dak_dvi = sympy.diff(ak, v_i)
# dak_dvj = sympy.diff(ak, v_j)
# dak_dvk = sympy.diff(ak, v_k)


# dai_dpi_latex = sympy.latex(dai_dpi)
# print(dai_dpi_latex)


"""
0.5*A*Cd*rho_0*w*(-p_i*w + v_j)*(p_j*w + v_i)*exp((R_E + h0 - r/H)/(m*v_rel) + 0.5*A*Cd*p_i*rho_0*(p_j*w + v_i)*v_rel*exp((R_E + h0 - r)/H)/(H*m*r)



"""
