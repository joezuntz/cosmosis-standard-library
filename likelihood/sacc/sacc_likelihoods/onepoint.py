import numpy as np
from cosmosis.datablock import BlockError
from scipy.interpolate import interp1d
import pathlib
import sys

def extract_one_point_prediction(sacc_data, block, data_type, section, **kwargs):
    category = kwargs.get("category")

    tracer_tuples = sacc_data.get_tracer_combinations(data_type)
    theory_vector = []
    observable_min_vector = []
    observable_max_vector= []
    bins_vector = []

    for t in tracer_tuples:
        assert len(t) == 1, "One-point likelihoods only support single tracer data types"
        
        # This should be an n(z) tracer
        tracer = t[0]
        b = int(tracer.split("_")[1]) + 1
        
        # This selection below needs to be generalised more 1pt statistics we implement (cluster richness for instance, ...).
        if category == "one_point_mass":
            x_theory = block[section, f"mass_{b}"]
        else:
            x_theory = block[section, f"luminosity_{b}"]
        theory = block[section, f"bin_{b}"]
        theory_spline = interp1d(x_theory, theory, bounds_error=False, fill_value="extrapolate")

        window = None
        for w in sacc_data.get_tag("windows", data_type, t):
            if (window is not None) and (w is not window):
                raise ValueError("Sacc likelihood currently assumes data types share a window object")
            window = w

        # TO-DO: Check if the window thing is ok for 1pt stats and how to do the binning here either way.
        if category == "one_point_mass":
            x_nominal = np.array(sacc_data.get_tag("mass", data_type, t))
        else:
            x_nominal = np.array(sacc_data.get_tag("luminosity", data_type, t))
        
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
        if category == "one_point_mass":
            observable_min_vector.append(sacc_data.get_tag("mass_min", data_type, t))
            observable_max_vector.append(sacc_data.get_tag("mass_max", data_type, t))
        else:
            observable_min_vector.append(sacc_data.get_tag("luminosity_min", data_type, t))
            observable_max_vector.append(sacc_data.get_tag("luminosity_max", data_type, t))

        bins_vector.append(np.repeat(b, len(binned_theory)))

    theory_vector = np.concatenate(theory_vector)

    if category == "one_point_mass":
        metadata_kwargs = {
            "min": "mass_min",
            "max": "mass_max",
            "bins": "mass_bin",
        }
    else:
        metadata_kwargs = {
            "min": "luminosity_min",
            "max": "luminosity_max",
            "bins": "luminosity_bin",
        }
    metadata = {
        metadata_kwargs["min"]: np.concatenate(observable_min_vector),
        metadata_kwargs["max"]: np.concatenate(observable_max_vector),
        metadata_kwargs["bins"]: np.concatenate(bins_vector),
    }
    return theory_vector, metadata





    


