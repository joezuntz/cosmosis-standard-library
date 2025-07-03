import numpy as np
from cosmosis.datablock import BlockError
import pathlib
import sys

# Get the SpectrumInterp class from the spec_tools module.
# Should really put this somewhere else!
twopoint_dir = pathlib.Path(__file__).parent.parent.parent.resolve() / "2pt"
print(twopoint_dir)
sys.path.append(str(twopoint_dir))
from spec_tools import SpectrumInterp


def extract_one_point_prediction(sacc_data, block, data_type, section, **kwargs):
    # data_type = galaxy_stellarmassfunction_hist
    # section = smf
    category = kwargs.get("category")

    tracer_tuples = sacc_data.get_tracer_combinations(data_type)
    theory_vector = []
    mass_min_vector = []
    mass_max_vector= []
    bins_vector = []
    for t in tracer_tuples:
        assert len(t) == 1, "One-point likelihoods only support single tracer data types"
        
        # This should be an n(z) tracer
        tracer = t[0]
        b = int(tracer.split("_")[1]) + 1
        
        if category == "one_point_mass":
            x_theory = block[section, f"mass_{b}"]
        else:
            x_theory = block[section, f"lum_{b}"]
        theory = block[section, f"bin_{b}"]
        theory_spline = SpectrumInterp(x_theory, theory)

        window = None
        w = sacc_data.get_tag("window", data_type, t)
        if (window is not None) and (w is not window):
            raise ValueError("Sacc likelihood currently assumes data types share a window object")
            window = w

        # TO-DO: Check if the window thing is ok for 1pt stats and how to do the binning here either way.
        if category == "one_point_mass":
            x_nominal = np.array(sacc_data.get_tag("mass", data_type, t))
        else:
            x_nominal = np.array(sacc_data.get_tag("lum", data_type, t))
        if window is None:
            binned_theory = theory_spline(x_nominal)
        else:
            x_window = window.values
            theory_interpolated = theory_spline(x_window)
            index = sacc_data.get_tag("window_ind", data_type, t)
            weight = window.weight[:, index]

            # The weight away should hopefully sum to 1 anyway but we should
            # probably not rely on that always being true.
            # TO-DO: Check this for real statistics, but should be ok.
            binned_theory = (weight @ theory_interpolated) / weight.sum()

        theory_vector.append(binned_theory)

        mass_min_vector.append(sacc_data.get_tag(f"mass_min_{b}", data_type, t))
        mass_max_vector.append(sacc_data.get_tag(f"mass_max_{b}", data_type, t))
        bins_vector.append(np.repeat(b, len(binned_theory)))

    theory_vector = np.concatenate(theory_vector)

    metadata = {
        #"mass_min": np.concatenate(mass_min_vector),
        #"mass_max": np.concatenate(mass_max_vector),
        "mass_bin": np.concatenate(bins_vector),
    }

    return theory_vector, metadata





    


