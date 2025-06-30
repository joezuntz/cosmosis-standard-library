import numpy as np

def extract_one_point_prediction(sacc_data, block, data_type, section):
    # data_type = galaxy_stellarmassfunction_hist_log
    # section = smf
    tracer_tuples = sacc_data.get_tracer_combinations(data_type)
    theory = []
    mass_min = []
    mass_max = []
    bins = []
    for t in tracer_tuples:
        assert len(t) == 1, "One-point likelihoods only support single tracer data types"
        
        # This should be an n(z) tracer
        tracer = t[0]
        b = int(tracer.split("_")[1])
        
        theory_log_mass = block[section, "log_mass"]
        theory_smf = block[section, f"bin_{b}"]

        obs_mass = sacc_data.get_tag("mass", data_type, t)
        obs_mass = np.array(obs_mass)

        theory_mass = 10 ** theory_log_mass
        smf_interp = np.interp(obs_mass, theory_mass, theory_smf)

        theory.append(smf_interp)

        mass_min.append(sacc_data.get_tag("min_mass", data_type, t))
        mass_max.append(sacc_data.get_tag("max_mass", data_type, t))
        bins.append(np.repeat(b, len(obs_mass)))


    metadata = {
        "mass_min": np.concatenate(mass_min),
        "mass_max": np.concatenate(mass_max),
        "mass_bin": np.concatenate(bins),
    }
    theory_vector = np.concatenate(theory)

    return theory_vector, metadata





    


