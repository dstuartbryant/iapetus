"""Error plotting module."""

from typing import List

from plotly import graph_objects as go
from plotly.subplots import make_subplots
from pydantic import BaseModel, root_validator, validator

SUBPLOT_TYPES = ["vertical", "horizontal", "grid"]


class ErrorPlottingError(Exception):
    pass


class ErrorStdInput(BaseModel):
    x: List[float]
    y_error: List[float]
    y_std: List[float]
    x_axis_title: str
    y_axis_title: str
    subplot_index: int = None
    subplot_type: str = None

    @validator("subplot_type")
    def check_subplot_type_option(cls, v):
        if not isinstance(v, type(None)):
            if v not in SUBPLOT_TYPES:
                raise ErrorPlottingError(
                    f"Unexpected `subplot_type` found: {v}"
                )
        return v

    @root_validator
    def both_subplot_options(cls, values):
        if (
            "subplot_index" in values.keys()
            and "subplot_type" not in values.keys()
        ):
            raise ErrorPlottingError(
                "`subplot_type` must be provided if `subplot_index` is used."
            )
        return values

    def _std(self, mode: str):
        if mode not in ["upper", "lower"]:
            raise ErrorPlottingError(
                f"Unexpected standard deviation mode found: {mode}"
            )
        if mode == "upper":
            coeff = 1
        else:
            coeff = -1
        return [
            self.y_error[idx] + coeff * self.y_std[idx]
            for idx in range(len(self.y_error))
        ]

    @property
    def x_axis_num(self):
        if self.subplot_type == "vertical":
            return self.subplot_index - 1

    @property
    def y_axis_num(self):
        if self.subplot_type == "vertical":
            return self.subplot_index - 1

    @property
    def x_axis_key(self):
        axis_num = self.x_axis_num
        key = "xaxis"
        if axis_num > 0:
            key += str(axis_num)
        return key

    @property
    def y_axis_key(self):
        axis_num = self.y_axis_num
        key = "yaxis"
        if axis_num > 0:
            key += str(axis_num)
        return key

    @property
    def upper_std(self):
        return self._std("upper")

    @property
    def lower_std(self):
        return self._std("lower")


class Plot1dErrorWithStd:
    """Class for plotting 1-D errors with their standard deviations."""

    def __init__(
        self,
        upper: bool,
        lower: bool,
    ):
        """
        Args:
            upper (bool): If True, plots standard deviations as upper bound of
                error
            lower (bool): If True, plots standard deviations as lower bound of
                error
        """
        self.upper = upper
        self.lower = lower
        self.figure = go.Figure()

    def add_error_trace(self, x, y):
        self.figure.add_trace(
            go.Scatter(
                x=x,
                y=y,
                line={"color": "black", "width": 1, "dash": "solid"},
                showlegend=False,
                name="Position error [m]",
            ),
        )

    def add_std_trace(self, x, y):
        self.figure.add_trace(
            go.Scatter(
                x=x,
                y=y,
                line={"color": "red", "width": 1, "dash": "dash"},
                showlegend=False,
                name=r"$1\sigma \text{ [m]}$",
            )
        )

    def __call__(self, data: ErrorStdInput):
        self.add_error_trace(data.x, data.y_error)
        if self.upper:
            self.add_std_trace(data.x, data.upper_std)
        if self.lower:
            self.add_std_trace(data.x, data.lower_std)

        self.figure.update_yaxes(title_text=data.y_axis_title)
        self.figure.update_xaxes(title_text=data.x_axis_title)
        self.figure.show()


class Plot3dErrorWithStdInput(BaseModel):
    """Model for call input of `Plot3dErrorWithStd` class."""

    top: ErrorStdInput
    middle: ErrorStdInput
    bottom: ErrorStdInput

    def __init__(self, **data):
        super().__init__(**data)
        self.top.subplot_index = 1
        self.top.subplot_type = "vertical"
        self.middle.subplot_index = 2
        self.middle.subplot_type = "vertical"
        self.bottom.subplot_index = 3
        self.bottom.subplot_type = "vertical"


class Plot3dErrorWithStd:
    """Class for plotting 3-D errors with their standard deviations."""

    def __init__(
        self,
        upper: bool,
        lower: bool,
        horizontal_spacing: float = 0.08,
        vertical_spacing: float = 0.08,
    ):
        """
        Args:
            upper (bool): If True, plots standard deviations as upper bound of
                error
            lower (bool): If True, plots standard deviations as lower bound of
                error
        """
        self.upper = upper
        self.lower = lower
        self.std_legend_count = 0
        self.figure = make_subplots(
            rows=3,
            cols=1,
            shared_xaxes=True,
            shared_yaxes=False,
            horizontal_spacing=horizontal_spacing,
            vertical_spacing=vertical_spacing,
        )

    def add_error_trace(self, x, y, plot_idx):
        showlegend = False
        if plot_idx == 1:
            showlegend = True
        self.figure.add_trace(
            go.Scatter(
                x=x,
                y=y,
                line={"color": "black", "width": 1, "dash": "solid"},
                showlegend=showlegend,
                name="Position error [m]",
            ),
            row=plot_idx,
            col=1,
        )

    def add_std_trace(self, x, y, plot_idx):
        showlegend = False
        if plot_idx == 1 and self.std_legend_count == 0:
            showlegend = True
            self.std_legend_count += 1
        self.figure.add_trace(
            go.Scatter(
                x=x,
                y=y,
                line={"color": "red", "width": 1, "dash": "dash"},
                showlegend=showlegend,
                name=r"$1\sigma \text{ [m]}$",
            ),
            row=plot_idx,
            col=1,
        )

    def __call__(self, data: Plot3dErrorWithStdInput):
        params = [("top", 1), ("middle", 2), ("bottom", 3)]
        for param in params:
            row_data = getattr(data, param[0])

            self.add_error_trace(
                getattr(row_data, "x"),
                getattr(row_data, "y_error"),
                plot_idx=param[1],
            )
            if self.upper:
                self.add_std_trace(
                    getattr(row_data, "x"),
                    getattr(row_data, "upper_std"),
                    plot_idx=param[1],
                )
            if self.lower:
                self.add_std_trace(
                    getattr(row_data, "x"),
                    getattr(row_data, "lower_std"),
                    plot_idx=param[1],
                )

            y_axis_title = getattr(row_data, "y_axis_title")
            self.figure.update_yaxes(
                title_text=y_axis_title, row=param[1], col=1
            )

        x_axis_title = getattr(row_data, "x_axis_title")
        self.figure.update_xaxes(title_text=x_axis_title, row=param[1], col=1)
        self.figure.show()

        return self.figure
