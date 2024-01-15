"""Two-Body with atmospheric drag propagation demo module."""

import json
from os import path

import numpy as np
from plotly import graph_objects as go
from plotly.subplots import make_subplots
from pydantic import BaseModel

from iapetus.propagation.dynamics.nonlinear.astro.xforms import ntw_matrix
from iapetus.propagation.propagators import AstroInit, AstroProp

CURR_DIR = path.dirname(path.abspath(__file__))


def add_error_component_trace(fig, t, error, row):
    fig.add_trace(
        go.Scatter(
            x=t,
            y=error,
            line={"color": "black", "width": 1, "dash": "solid"},
            showlegend=False,
        ),
        row=row,
        col=1,
    )
    return fig


class StateVector(BaseModel):
    position_i_m: float
    position_j_m: float
    position_k_m: float
    velocity_i_mps: float
    velocity_j_mps: float
    velocity_k_mps: float


class State(BaseModel):
    timestamp: float
    vector: StateVector


fpath = path.join(CURR_DIR, "one_orbit_data.json")
data = json.load(open(fpath, "r"))

# T = [x["timestamp"] for x in data]
T = np.arange(0, 2 * 90 * 60)
first_state = State(**data[0])

t0 = T[0]
X0_dict = first_state.vector.dict()

# --------------------- Propagate without drag -----------------
tspan = T
dt = T[1] - T[0]
tspantol = 1e-1

a_init = AstroInit(
    state_vector_content=["translational"],
    celestial_body="Earth",
    stm_flag=False,
    integrator="rk45",
    perturbations=[],
)
aprop = AstroProp(a_init)

T, Y_no_drag, _ = aprop(tspan, dt, tspantol, ui_state=X0_dict)

# --------------------- Propagate WITH drag -----------------

X0_dict["A_m2"] = 3**2
X0_dict["Cd"] = 2.0
X0_dict["m_kg"] = 50

a_drag_init = AstroInit(
    state_vector_content=["translational"],
    celestial_body="Earth",
    stm_flag=False,
    integrator="rk45",
    perturbations=["atmospheric-drag"],
)

aprop_drag = AstroProp(a_drag_init)

T, Y_drag, _ = aprop_drag(tspan, dt, tspantol, ui_state=X0_dict)

# ---------------------- Compute RIC Errors ---------------
e_r = []
e_i = []
e_c = []

for idx, y_no_drag in enumerate(Y_no_drag):
    y_drag = Y_drag[idx]

    p_error = y_no_drag[:3] - y_drag[:3]

    R_ECI_to_NTW = ntw_matrix(y_no_drag[3:], y_no_drag[:3])
    p_error_ntw = R_ECI_to_NTW @ p_error

    e_r.append(p_error_ntw[0])
    e_i.append(p_error_ntw[1])
    e_c.append(p_error_ntw[2])

# ----------------------- Plot --------------------------

fig = make_subplots(
    rows=3,
    cols=1,
    shared_xaxes=True,
    shared_yaxes=False,
    horizontal_spacing=0.08,
    vertical_spacing=0.08,
)
fig = add_error_component_trace(fig, T, e_r, 1)
fig["layout"]["yaxis"]["title"] = "Radial error [m]"

fig = add_error_component_trace(fig, T, e_i, 2)
fig["layout"]["yaxis2"]["title"] = "In-Track error [m]"

fig = add_error_component_trace(fig, T, e_c, 3)
fig["layout"]["yaxis3"]["title"] = "Cross-Track error [m]"

fig["layout"]["xaxis3"]["title"] = "Time elapsed past epoch [s]"
fig.show()
