from builtins import range
def setup(options):
    return {}


def execute(block, config):
    bias_section = "bin_bias"

    section = "galaxy_cl"
    nbin = block[section, "nbin"]
    for i in range(nbin):
        b1 = block[bias_section, "b_{}".format(i + 1)]
        for j in range(nbin):
            b2 = block[bias_section, "b_{}".format(i + 1)]
            name = "bin_{}_{}".format(i + 1, j + 1)
            if block.has_value(section, name):
                block[section, name] *= b1 * b2
    return 0
