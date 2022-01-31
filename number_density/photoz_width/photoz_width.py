from builtins import range
from cosmosis.datablock import option_section, names
from scipy.interpolate import interp1d
import numpy as np

MODES = ["stretch"]


def setup(options):
    mode = options.get_string(option_section, "mode", "stretch")
    sample = options.get_string(option_section, "sample", "")
    interpolation = options.get_string(
        option_section, "interpolation", "cubic")
    bias_section = options.get_string(option_section, "bias_section", "")
    per_bin = options.get_bool(option_section, "per_bin", True)
    if sample == "":
        pz = names.wl_number_density
    else:
        pz = sample
    if bias_section == "" and sample == "":
        bias_section = "wl_photoz_errors"
    elif bias_section == "":
        bias_section = sample + "_errors"
    if mode not in MODES:
        raise ValueError("mode for photoz must be one of: %r" % MODES)
    return {"mode": mode, "sample": pz, "bias_section": bias_section, "interpolation": interpolation, "per_bin": per_bin}


def execute(block, config):
    mode = config['mode']
    pz = config['sample']
    interpolation = config['interpolation']
    biases = config['bias_section']
    nbin = block[pz, "nbin"]
    z = block[pz, "z"]
    for i in range(1, nbin + 1):
        bin_name = "bin_%d" % i
        nz = block[pz, bin_name]
        if config["per_bin"]:
            stretch = block[biases, "width_%d" % i]
        else:
            stretch = block[biases, "width_0"]

        if mode == "stretch":
            zmean = np.average(z,weights=nz)
            f = interp1d(stretch*(z-zmean)+zmean, nz, kind=interpolation, 
                fill_value=0.0, bounds_error=False)
            nz_biased = f(z)
        else:
            raise ValueError("Unknown photo-z mode")

        # normalize
        nz_biased /= np.trapz(nz_biased, z)
        block[pz, bin_name] = nz_biased
    return 0


def cleanup(config):
    pass
