from cosmosis.datablock import option_section, names

def setup(options):
    suffix = options.get_string(option_section, "suffix", "")
    return suffix

def execute(block, config):
    suffix = config

    for sec,name in block.keys('intrinsic_alignment_parameters'):
        if (sec=='intrinsic_alignment_parameters') and (suffix in name):
            out_name = name.replace('_%s'%suffix,'')
            block[sec,out_name] = block[sec,name]

    return 0
