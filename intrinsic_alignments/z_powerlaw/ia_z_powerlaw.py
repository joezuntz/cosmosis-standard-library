from cosmosis.datablock import names, option_section
import numpy as np

def setup(options):
    name = options.get_string(option_section, "name", default="").lower()
    if name:
        suffix = "_" + name
    else:
        suffix = ""
    return {"suffix":suffix}

def execute(block, config):
    #Get the names of the sections to save to
    suffix = config['suffix']
    ia_section = names.intrinsic_alignment_parameters
    ia_ii = names.intrinsic_power + suffix
    ia_mi = names.matter_intrinsic_power + suffix

    #read in power spectra for P_II and P_MI
    z,k,p_ii=block.get_grid(ia_ii,"z","k_h","p_k")
    z,k,p_mi=block.get_grid(ia_mi,"z","k_h","p_k")

    #read alpha from ia_section values section
    alpha = block[ia_section,'alpha']
    _,z_grid=np.meshgrid(k,z)

    #Construct and apply redshift scaling
    z_scaling=(1+z_grid)**alpha
    p_ii*=z_scaling**2
    p_mi*=z_scaling

    #Save grid back to the block
    block.replace_grid(ia_ii, "z", z, "k_h", k, "p_k", p_ii)
    block.replace_grid(ia_mi, "z", z, "k_h", k, "p_k", p_mi)

    return 0

def cleanup(config):
    pass
