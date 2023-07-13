# Documentation on how to use derivedParameters.py can be found on the SFB
# confluence
def computeDeltaT(values):
    import math

    # Give explicitly prescribed time step sizes precedence over
    # computed values
    if "delta_t" in values:
        return values["delta_t"]

    domainLength = 3e-3
    resolution = values["resolution"]

    h = domainLength/(resolution)

    # Check for the CFL number
    cfl_factor = values["CFL_num"]
    u_characteristic = 30
    cfl_delta_t = cfl_factor*h/u_characteristic

    return cfl_delta_t

delta_t = computeDeltaT(locals())

# This dummy is required to avoid a PyFoam error if delta_t is given explicitly in the
# parameter file
dummy = 0
