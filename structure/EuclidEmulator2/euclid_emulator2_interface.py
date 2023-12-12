import sys
import os
from cosmosis.datablock import names, option_section


dirname = os.path.split(__file__)[0]

# enable debugging from the same directory
if not dirname.strip():
    dirname = '.'



def setup(options):
    input_section = options.get_string(option_section, "input_section", default=names.matter_power_lin)
    output_section = options.get_string(option_section, "output_section", default=names.matter_power_nl)

    # Set up the path to let us import the emulator
    install_dir = f"{dirname}/ee_install/"
    sys.path.insert(0, install_dir)

    # check everything imports
    import euclidemu2
    print("Loaded Euclid Emulator 2 installation at", euclidemu2.__file__)
    return [input_section, output_section]


def execute(block, config):
    import euclidemu2 as ee2

    # Recover config information
    input_section, output_section = config

    # Get cosmo params from block
    pars = names.cosmological_parameters
    params = {
        'As': block[pars, "A_s"],
        'ns': block[pars, "n_s"],
        'Omb': block[pars, "Omega_b"],
        'Omm': block[pars, "Omega_m"],
        'h': block[pars, "h0"],
        'mnu': block[pars, "mnu"],
        'w': block[pars, "w"],
        'wa': block[pars, "wa"],
    }

    # Get z and k from the NL power section
    z, k, P = block.get_grid(input_section, "z", "k_h", "P_k")

    if len(z) > 100:
        raise ValueError("EuclidEmulator2 only allows up to 100 redshift values")

    _, b = ee2.get_boost(params, z, k)

    # Not sure why but b comes back as a dictionary of arrays
    # instead of a 2D array
    P_boosted = P.copy()
    for i in range(len(z)):
        P_boosted[i] *= b[i]

    if output_section == input_section:
        # save the original grid to a new section, with _dm added
        # and save the new grid back to the input
        block.put_grid(input_section + "_dm", "z", z, "k_h", k, "P_k", P)
        block.replace_grid(input_section, "z", z, "k_h", k, "P_k", P_boosted)
    else:
        # just save the output
        block.put_grid(output_section, "z", z, "k_h", k, "P_k", P_boosted)

    return 0

if __name__ == '__main__':
    setup(None)