from __future__ import print_function
from builtins import range
from cosmosis.datablock import option_section, names


def setup(options):
    do_shear_shear = options.get_bool(option_section, "shear-shear", True)
    do_position_shear = options.get_bool(
        option_section, "position-shear", True)
    perbin = options.get_bool(option_section, "perbin", False)

    suffix = options.get_string(option_section, "suffix", "")

    print()
    print("The add_intrinsic module will try to combine IA terms into")
    if do_shear_shear and do_position_shear:
        print("both the shear-shear and position-shear spectra.")
    elif do_shear_shear:
        print("only the shear-shear spectra.")
    elif do_position_shear:
        print("only the position-shear.")
    else:
        print("... actually not into anything. You set shear-shear=F and position-shear=F")
        print("Ths module will not do anything in this configuration")
    print()

    if suffix:
        suffix = "_" + suffix

    sec_names = {
        "shear_shear": "shear_cl" + suffix,
        "shear_shear_bb": "shear_cl_bb" + suffix,
        "shear_shear_gg": "shear_cl_gg" + suffix,
        "galaxy_shear": "galaxy_shear_cl" + suffix,
        "galaxy_intrinsic": "galaxy_intrinsic_cl"  + suffix,
        "intrinsic_intrinsic": "shear_cl_ii"  + suffix,
        "intrinsic_intrinsic_bb": "shear_cl_ii_bb"  + suffix,
        "parameters": "intrinsic_alignment_parameters" + suffix,
    }   

    return do_shear_shear, do_position_shear, perbin, sec_names


def execute(block, config):
    do_shear_shear, do_position_shear, perbin, sec_names = config

    shear_shear = sec_names['shear_shear']
    shear_shear_bb = sec_names['shear_shear_bb']
    shear_shear_gg = sec_names['shear_shear']
    galaxy_shear = sec_names['galaxy_shear']
    galaxy_intrinsic = sec_names['galaxy_intrinsic']
    parameters = sec_names['parameters']
    intrinsic_intrinsic = sec_names['intrinsic_intrinsic']
    intrinsic_intrinsic_bb = sec_names['intrinsic_intrinsic_bb']

    if do_shear_shear:
        nbin_shear = block[shear_shear, 'nbin']
    elif do_position_shear:
        nbin_shear = block[galaxy_intrinsic, 'nbin_b']
    if do_position_shear:
        nbin_pos = block[galaxy_shear, 'nbin_a']

    if perbin:
        A = [block[parameters, "A{}".format(i + 1)]
             for i in range(nbin_shear)]
    else:
        A = [1 for i in range(nbin_shear)]

    if do_shear_shear:
        # for shear-shear, we're replacing 'shear_cl' (the GG term) with GG+GI+II...
        # so in case useful, save the GG term to shear_cl_gg.
        # also check for a b-mode contribution from IAs
        block[shear_shear_gg, 'ell'] = block[shear_shear, 'ell']
        for i in range(nbin_shear):
            for j in range(i + 1):
                bin_ij = 'bin_{0}_{1}'.format(i + 1, j + 1)
                bin_ji = 'bin_{1}_{0}'.format(i + 1, j + 1)
                block[shear_shear_gg, bin_ij] = block[shear_shear, bin_ij]
                block[shear_shear, bin_ij] += (
                    A[i] * A[j] * block[intrinsic_intrinsic, bin_ij]  # II
                    + A[j] * block[galaxy_intrinsic,
                                   bin_ij]  # The two GI terms
                    + A[i] * block[galaxy_intrinsic, bin_ji]
                )
                if block.has_section(intrinsic_intrinsic_bb):
                    block[shear_shear_bb, bin_ij] = block[intrinsic_intrinsic_bb, bin_ij]
    if do_position_shear:
        for i in range(nbin_pos):
            for j in range(nbin_shear):
                bin_ij = 'bin_{0}_{1}'.format(i + 1, j + 1)
                block[galaxy_shear, bin_ij] += A[j] * \
                    block[galaxy_intrinsic, bin_ij]
    return 0
