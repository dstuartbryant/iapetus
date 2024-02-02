"""Integrators module."""

from pydantic import BaseModel, Field


class IntegrateConfig(BaseModel):
    """Integrator configuration class."""

    dt: float = Field(..., description="Integrator time step size")
    tol: float = Field(
        ...,
        description=(
            "integrator time step tolerance; used to verify that time steps "
            "that divide an integration time span are divided evenly. "
            "Typically, a good value is 1e-6, but if numerical issues are "
            "experienced, or, if faster compute times are desired at the cost "
            "of accuracy, try a larger value."
        ),
    )
