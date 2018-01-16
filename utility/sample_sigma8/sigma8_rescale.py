from numpy import log, pi
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section

cosmo = section_names.cosmological_parameters
cmb = section_names.cmb_cl
matter_powspec = section_names.matter_power_lin


def setup(options):
    return 0


def execute(block, config):

    # Get parameters from sampler and CAMB output
    sigma8_input = block[cosmo, 'sigma8_input']
    sigma8_camb = block[cosmo, 'sigma_8']

    A_s = block[cosmo, 'A_s']
    TT = block[cmb, 'TT']
    EE = block[cmb, 'EE']
    BB = block[cmb, 'BB']
    TE = block[cmb, 'TE']
    P_k = block[matter_powspec, 'P_k']

    zmin = block[matter_powspec, 'z'].min()
    if zmin != 0.0:
        raise ValueError(
            "You need to set zmin=0 in CAMB to use the sigma8_rescale module.")

    # Calculate rescale factor
    r = (sigma8_input**2) / (sigma8_camb**2)

    # Rescale CMB Cl and matter power outputs
    A_s *= r
    TT *= r
    EE *= r
    BB *= r
    TE *= r
    P_k *= r

    # Save back to block
    block[cosmo, 'A_s'] = A_s
    block[cmb, 'TT'] = TT
    block[cmb, 'EE'] = EE
    block[cmb, 'BB'] = BB
    block[cmb, 'TE'] = TE
    block[matter_powspec, 'P_k'] = P_k

    block[cosmo, 'sigma_8'] = sigma8_input

    # signal that everything went fine
    return 0


def cleanup(config):
    # nothing to do here!  We just include this
    # for completeness
    return 0
