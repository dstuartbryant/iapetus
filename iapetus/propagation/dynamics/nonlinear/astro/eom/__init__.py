"""Astrodynamics equations of motion (EOM) module."""


# class Eom:
#     """Astrodynamics equations of motion (EOM) class."""

#     def __init__(
#         self,
#         param: EomParam,
#         eom_list: List[str] = None,
#     ):
#         """
#         Args:
#             param (dict): Parameters dictionary
#             eom_list (List[str]): List of EOM dynamics to use, defaults to None
#         """
#         self.param = param
#         if eom_list:
#             self.eom_list = eom_list
#         else:
#             self.eom_list = []
#         self._validate_params()
#         self.dynamics = self._configure()
