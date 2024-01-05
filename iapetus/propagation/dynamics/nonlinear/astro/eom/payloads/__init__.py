"""Astrodynamics equations of motion (EOM) payloads module.

'Payload' refers to objects that travel 'through' methods, where the methods
call specific attributes as needed. This is an attempt at unifiying interfaces
to minimize code.

"""


from .configs import AtmosphericDragInitConfig, TwoBodyInitConfig
