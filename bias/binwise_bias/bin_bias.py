from cosmosis.datablock import names, option_section
import sys

def setup(options):
    perbin=options.get_bool("bin_bias","perbin",True)
    auto_only=options.get_bool("bin_bias","auto_only",False)
    return perbin,auto_only

def execute(block, options):
    perbin,auto_only=options
    if block.has_section('galaxy_shear_cl'):
        n_z_bins_shear=block["galaxy_shear_cl","nbin_b"]
        n_z_bins_pos=block["galaxy_shear_cl","nbin_a"]
        apply_to_cl=True
    elif block.has_section('galaxy_cl'):
        n_z_bins_pos=block["galaxy_cl","nbin"]
        apply_to_cl=True
    elif block.has_section('galaxy_shear_xi'):
        n_z_bins_shear=block["galaxy_shear_xi","nbin_b"]
        n_z_bins_pos=block["galaxy_shear_xi","nbin_a"]
        apply_to_cl=False
    elif block.has_section('galaxy_xi'):
        n_z_bins_pos=block["galaxy_xi","nbin"]
        apply_to_cl=False
    else:
        sys.stderr.write("ERROR: The bin_bias module could not find any of galaxy_cl, galaxy_shear_cl, galaxy_shear_xi, or galaxy_xi to bias\n")
        return 1

    #We may be doing per-bin biases or a single global value
    if perbin:
        #per-bin - use b0,b1,b2, ...
        biases = [block["bin_bias", "b%d"%pos_bin] for pos_bin in xrange(1,n_z_bins_pos+1)]
    else:
        #all the same - just use b0
        biases = [block["bin_bias", "b0"] for pos_bin in xrange(n_z_bins_pos)]

    for pos_bin1 in xrange(n_z_bins_pos):
        bias1 = biases[pos_bin1]

        if apply_to_cl:
            if block.has_section('galaxy_shear_cl'):
                for shear_bin in xrange(n_z_bins_shear):
                    name = "bin_{}_{}".format(pos_bin1+1, shear_bin+1)
                    block["galaxy_shear_cl", name] *= bias1

            if block.has_section('galaxy_intrinsic_cl'):
                for shear_bin in xrange(n_z_bins_shear):
                    name = "bin_{}_{}".format(pos_bin1+1, shear_bin+1)
                    block["galaxy_intrinsic_cl", name] *= bias1

            if block.has_section('galaxy_cl'):
                for pos_bin2 in xrange(pos_bin1+1):
                    bias2 = biases[pos_bin2]
                    name = "bin_{}_{}".format(pos_bin1+1, pos_bin2+1)    
                    block['galaxy_cl', name] *= bias1*bias2

        else:
            if block.has_section('galaxy_xi'):
                for pos_bin2 in xrange(pos_bin1+1):
                    if auto_only:
                        if pos_bin2!=pos_bin1:
                            continue
                    bias2 = biases[pos_bin2]
                    name = "bin_{}_{}".format(pos_bin1+1, pos_bin2+1)    
                    block['galaxy_xi', name] *= bias1*bias2            

            if block.has_section('galaxy_shear_xi'):
                for shear_bin in zrange(n_z_bins_shear):
                    name = "bin_{}_{}".format(pos_bin1+1, shear_bin+1)    
                    block['galaxy_shear_xi', name] *= bias1            


    return 0
