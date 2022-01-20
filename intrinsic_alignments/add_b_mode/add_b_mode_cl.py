from cosmosis.datablock import option_section, names

def setup(options):
    ee_section=options.get_string(option_section, "ee_section", "shear_cl")
    bb_section=options.get_string(option_section, "bb_section", "shear_cl_bb")
    return ee_section, bb_section

def execute(block, config):
    ee_section, bb_section = config

    nbin_shear = block[ee_section, 'nbin']
    p_section, m_section = "shear_cl_eplusb", "shear_cl_eminusb"

    # clone the EE section, such that we retain all of the metadata
    block._copy_section(ee_section,p_section)
    block._copy_section(ee_section,m_section)

    ell = block[ee_section, "ell"]
    block[p_section, "ell"] = ell
    block[m_section, "ell"] = ell
    block[p_section, "nbin"] = nbin_shear
    block[m_section, "nbin"] = nbin_shear
    for i in range(nbin_shear):
        for j in range(0,i+1):
            bin_ij = 'bin_{0}_{1}'.format(i+1,j+1)
            ee = block[ee_section, bin_ij]
            bb = block[bb_section, bin_ij]
            block[p_section, bin_ij] = ee+bb
            block[m_section, bin_ij] = ee-bb

    return 0
