# Documentation on how to use derivedParameters.py can be found on the SFB
# confluence
def computeDeltaT(values):
    import math

    # Give explicitly prescribed time step sizes precedence over
    # computed values
    if "delta_t" in values:
        return values["delta_t"]

    domainLength = 16.0*values["radius"]
    resolution = values["n_base"] # ["n_base"]->block # ["resolution"]->hex
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
    ltria = domainLength/(resolution)

    # Safety coefficient to stay below critical delta_t threshold
    sf = values["scale_delta_t"]

    return sf*math.sqrt((rho_droplet + rho_ambient)*math.pow(ltria,3.0)/(2.0*math.pi*sigma))

def computeBoxLength(values):
    """Computes the dimension of the cubic bounding box for the Duineveld
    rising buble simulation in an ALE Relative Reference Frame (RRF) from the
    bubble radius.""" 

    radius = values["radius"]

    return 20*radius;


delta_t = computeDeltaT(locals())
box_length = computeBoxLength(locals())

