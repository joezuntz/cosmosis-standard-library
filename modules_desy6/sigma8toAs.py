################
#THIS MODULE IS A COPY FROM THE KIDS KCAP repository
#Tilman Troester committed the module
#Need to be tested and chceked for our purpose
################


from cosmosis.datablock import option_section, names

import camb
import numpy as np

def setup(options):
    sigma8_name = options.get_string(option_section, "sigma8_name", default="sigma_8_input")
    kmax = options.get_double(option_section, "kmax", default=1.2)
    k_per_logint = options.get_int(option_section, "k_per_logint", default=5)
    return sigma8_name, kmax, k_per_logint

def execute(block, config):
    sigma8_name, kmax, k_per_logint = config

    h = block[names.cosmological_parameters, "h0"]
    H0 = h*100
    ombh2 = block[names.cosmological_parameters, "ombh2"]
    omch2 = block[names.cosmological_parameters, "omch2"]
    omk = block[names.cosmological_parameters, "omega_k"]
    mnu = block.get_double(names.cosmological_parameters, "mnu", default=0.06)

    w = block.get_double(names.cosmological_parameters, "w", default=-1.0)
    wa = block.get_double(names.cosmological_parameters, "wa", default=0.0)

    ns = block[names.cosmological_parameters, "n_s"]
    sigma8 = block[names.cosmological_parameters, sigma8_name]

    fid_As = 2.1e-9

    p = camb.CAMBparams(WantTransfer=True,
                        Want_CMB=False, Want_CMB_lensing=False, DoLensing=False,
                        NonLinear="NonLinear_none",
                        WantTensors=False, WantVectors=False, WantCls=False,
                        WantDerivedParameters=False,
                        want_zdrag=False, want_zstar=False)
    p.set_accuracy(DoLateRadTruncation=True)
    p.Transfer.high_precision = False
    p.Transfer.accurate_massive_neutrino_transfers = False
    p.Transfer.kmax = kmax
    p.Transfer.k_per_logint = k_per_logint
    p.Transfer.PK_redshifts = np.array([0.0])

    p.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, omk=omk, mnu=mnu)
    p.set_dark_energy(w=w, wa=wa)
    p.set_initial_power(camb.initialpower.InitialPowerLaw(As=fid_As, ns=ns))

    p.Reion = camb.reionization.TanhReionization()
    p.Reion.Reionization = False

    r = camb.get_results(p)

    fid_sigma8 = r.get_sigma8()[-1]

    As = fid_As*(sigma8/fid_sigma8)**2

    print("fid_As",fid_As,'As',As,'sigma8',sigma8,'fid_sigma8',fid_sigma8)
    block[names.cosmological_parameters, "A_s"] = As

    return 0

def cleanup(config):
    pass
