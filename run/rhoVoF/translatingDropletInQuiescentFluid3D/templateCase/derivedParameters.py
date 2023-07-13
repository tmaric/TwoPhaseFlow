# Documentation on how to use derivedParameters.py can be found on the SFB
# confluence
from __future__ import division
def computeDeltaT(values):
    import math

    # Give explicitly prescribed time step sizes precedence over
    # computed values
    if "delta_t" in values:
        return values["delta_t"]

    domainLength = 8.0*values["radius"]
    resolution = values["resolution"]
    rho_droplet = values["rho_droplet"]

    # rho_ambient might be expressed in terms of rho_droplet
    rho_ambient = 0.0
    if values["rho_ambient"] == "$rho_droplet":
        rho_ambient = rho_droplet
    else:
        rho_ambient = values["rho_ambient"]
    sigma = values["surface_tension_coefficient"]

    # Use the resolution of the Eulerian mesh rather than the front element size.
    # So far there is no indication that the triangle size plays a role for the
    # maximum time step.
    # The computation follows eq. (18) in Popinet 2009 and
    # eq. in (43) Denner & van Wachem 2015, respectively.
    h = domainLength/(resolution)

    # Safety coefficient to stay below critical delta_t threshold
    sf = values["scale_delta_t"]
    if sigma != 0:
       capillary_delta_t = sf*math.sqrt((rho_droplet + rho_ambient)*math.pow(h,3.0)/(2.0*math.pi*sigma))
    else:
        capillary_delta_t = 10000
    # Check for the CFL number
    cfl_factor = values["cfl_factor"]
    u_characteristic = values["z_velocity"]
    cfl_delta_t = cfl_factor*h/u_characteristic

    return min(capillary_delta_t, cfl_delta_t)

delta_t = computeDeltaT(locals())

# This dummy is required to avoid a PyFoam error if delta_t is given explicitly in the
# parameter file
dummy = 0
