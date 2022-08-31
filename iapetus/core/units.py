"""Measurement unit tools."""

# Systems of measure
SYSSI = "SI"
SYSNONE = "NONE"


SYSTYPES = [SYSSI, SYSNONE]


"""
Don't like where this is going. Need to revisit this.

"""


# class MeasSystem:
#     "System of measurement base class."

#     def __init__(self, meas_system: str):
#         """
#         Args:
#             meas_system (str): indicates the system of measure.
#         """
#         if meas_system not in SYSTYPES:
#             raise ValueError(
#                 f"Unexpected measurement system type found: {meas_system}."
#             )
#         self.system = meas_system
#         self._configure()

#     def _configure(self):
#         if self.system == SYSNONE:
#             self.default_distance_unit = None
#             self.default_linear_speed_unit = None
#             self.default_time_unit = None


# class Unit:
#     """Unit base class."""

#     def __init__(self, unit_system: str):
#         """
#         Args:
#             unit_system (str): indicates the system used for specifying units
#                 of measure.
#         """
#         self.system = unit_system
