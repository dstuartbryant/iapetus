"""Two-Body propagation demo module."""

import json
from os import path

from plotly import graph_objects as go
from plotly.subplots import make_subplots
from pydantic import BaseModel

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

truth_data = json.load(open(fpath, "r"))
truth = [State(**x) for x in truth_data]

t0 = truth[0].timestamp
X0_dict = truth[0].vector.dict()

tspan = [x.timestamp for x in truth]
dt = tspan[1] - tspan[0]
tspantol = 1e-1

a_init = AstroInit(
    state_vector_content=["translational"],
    celestial_body="Earth",
    stm_flag=False,
    integrator="rk45",
    perturbations=[],
)
aprop = AstroProp(a_init)

T, Y, _ = aprop(tspan, dt, tspantol, ui_state=X0_dict)

e_i = []
e_j = []
e_k = []
for idx, y in enumerate(Y):
    e_i.append(truth[idx].vector.position_i_m - y[0])
    e_j.append(truth[idx].vector.position_j_m - y[1])
    e_k.append(truth[idx].vector.position_k_m - y[2])

fig = make_subplots(
    rows=3,
    cols=1,
    shared_xaxes=True,
    shared_yaxes=False,
    horizontal_spacing=0.08,
    vertical_spacing=0.08,
)
fig = add_error_component_trace(fig, T, e_i, 1)
fig["layout"]["yaxis"]["title"] = "I-axis error [m]"

fig = add_error_component_trace(fig, T, e_j, 2)
fig["layout"]["yaxis2"]["title"] = "J-axis error [m]"

fig = add_error_component_trace(fig, T, e_k, 3)
fig["layout"]["yaxis3"]["title"] = "K-axis error [m]"

fig["layout"]["xaxis3"]["title"] = "Time elapsed past epoch [s]"
fig.show()
