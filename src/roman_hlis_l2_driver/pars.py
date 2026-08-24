"""
General parameter information goes here.
"""

# PyIMCOM settings should be accessible from here.
from pyimcom.config import Settings as Stn

# global information about Roman
n_sca = 18

n_side = Stn.sca_nside

# if this isn't true, we'll have a problem
assert len(Stn.SCAFov) >= n_sca

# gain information
ref_gain = 1.458
