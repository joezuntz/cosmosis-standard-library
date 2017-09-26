from __future__ import print_function
from cosmosis.datablock import option_section, names


def setup(options):
    do_shear_shear = options.get_bool(option_section, "shear-shear", True)
    do_position_shear = options.get_bool(
        option_section, "position-shear", True)
    perbin = options.get_bool(option_section, "perbin", False)
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
    return do_shear_shear, do_position_shear, perbin


def execute(block, config):
    do_shear_shear, do_position_shear, perbin = config

    if do_shear_shear:
        nbin_shear = block[names.shear_cl, 'nbin']
    elif do_position_shear:
        nbin_shear = block["galaxy_intrinsic_cl", 'nbin_b']
    if do_position_shear:
        nbin_pos = block["galaxy_shear_cl", 'nbin_a']

    if perbin:
        A = [block[names.intrinsic_alignment_parameters, "A{}".format(i + 1)]
             for i in xrange(nbin_shear)]
    else:
        A = [1 for i in xrange(nbin_shear)]

    if do_shear_shear:
        # for shear-shear, we're replacing 'shear_cl' (the GG term) with GG+GI+II...
        # so in case useful, save the GG term to shear_cl_gg.
        # also check for a b-mode contribution from IAs
        block[names.shear_cl_gg, 'ell'] = block[names.shear_cl, 'ell']
        for i in xrange(nbin_shear):
            for j in xrange(i + 1):
                bin_ij = 'bin_{0}_{1}'.format(i + 1, j + 1)
                bin_ji = 'bin_{1}_{0}'.format(i + 1, j + 1)
                block[names.shear_cl_gg, bin_ij] = block[names.shear_cl, bin_ij]
                block[names.shear_cl, bin_ij] += (
                    A[i] * A[j] * block[names.shear_cl_ii, bin_ij]  # II
                    + A[j] * block[names.shear_cl_gi,
                                   bin_ij]  # The two GI terms
                    + A[i] * block[names.shear_cl_gi, bin_ji]
                )
                if block.has_section("shear_cl_ii_bb"):
                    block["shear_cl_bb", bin_ij] = block["shear_cl_ii_bb", bin_ij]
    if do_position_shear:
        for i in xrange(nbin_pos):
            for j in xrange(nbin_shear):
                bin_ij = 'bin_{0}_{1}'.format(i + 1, j + 1)
                block["galaxy_shear_cl", bin_ij] += A[j] * \
                    block["galaxy_intrinsic_cl", bin_ij]
    return 0
