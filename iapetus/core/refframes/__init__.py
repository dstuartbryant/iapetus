"""Reference frame handling module."""


class ReferenceFrame:
    """Reference frame base class.

    Attributes:
        name (str): uniquely identifies the reference frame.
    """

    def __init__(self, name: str):
        self.name = name


plane = ReferenceFrame("2D Plane")
