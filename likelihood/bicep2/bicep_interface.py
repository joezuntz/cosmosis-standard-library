from __future__ import print_function
import bicep_util as bu
from contextlib import contextmanager
import os
import sys
from cosmosis.datablock.cosmosis_py.section_names import likelihoods, cmb_cl

dirname = os.path.split(__file__)[0]

import numpy as np


@contextmanager
def working_directory(dirname):
    """This context manager is used with the 'with' statement to switch temporarily to another directory """
    # If the directory does not exist we first create it.
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    # Record the original starting directory so we can go back to it at the end
    start_dir = os.getcwd()
    # Move to the chosen new directory
    os.chdir(dirname)
    # Return to the calling code and run the statement in the "with" block.
    try:
        yield
    # Whether or not there is any kind of error in the "with" block we finally return to
    # the directory we started in.
    finally:
        os.chdir(start_dir)


def setup(options):
    # From the directory the data is stored in,
    # load the various files they want.
    print("Loading BICEP2 data from directory", dirname)

    with working_directory(dirname):
        # Call BICEP's setup function with the
        # list of what we want
        experiment = 'bicep2'
        field = 'B'
        lcdm = "B2_3yr_camb_planck_lensed_uK_20140314.txt"
        C_l, C_l_hat, N_l, C_fl, M_inv, bpwf_l, bpwf_Cs_l = bu.init(
            experiment, field)

        # BICEP needs two different fiducial power spectra, which we load from
        # file here
        (ell1, inpmodLCDM_Cs_l) = bu.load_cmbfast(lcdm)
        expvLCDM = bu.calc_expvals(ell1, inpmodLCDM_Cs_l, bpwf_l, bpwf_Cs_l)
        file_in = "B2_3yr_camb_planck_withB_uK_20140314.txt"
        (ell, BB_Cs_l) = bu.load_cmbfast(file_in)

        # A bit later on in this code in the execute function we implicitly
        # assume that ell starts at zero.  As this is the kind of thing
        # that might change we explicitly check for it here.
        assert int(
            ell[0]) == 0, "BICEP data format error - something has changed - see code"

        # A quick warning for the unwary user...
        print()
        print("You're using the BICEP2 likelihood code.")
        print("That's great, but just one quick note: ")
        print("If using camb for your spectra you need do_nonlinear = T, k_eta_max_scalar = 14000 to get")
        print("accurate enough lensing B-modes for good results.  There will be equivalents for other codes too.")
        print()
        print("If you're running a demo then don't worry, we took care of this for you.")
        print()
    return C_l, C_l_hat, N_l, C_fl, M_inv, bpwf_l, bpwf_Cs_l, ell1, inpmodLCDM_Cs_l, expvLCDM, BB_Cs_l, ell


def execute(block, config):
    # Get back all the stuff we loaded in during setup.
    C_l, C_l_hat, N_l, C_fl, M_inv, bpwf_l, bpwf_Cs_l, ell1, inpmodLCDM_Cs_l, expvLCDM, BB_Cs_l, ell = config

    # Load ell and Cls
    # from the package
    # BICEP code expects the following order: TT, TE, EE, BB, TB, EB, ET, BT, BE
    # Put theory Cls into a copy of what we loaded in from the init
    # to get the shape right
    this_mod = BB_Cs_l.copy()
    # Get right ell range
    ell_theory = block[cmb_cl, "ELL"]
    theory_ell_min = int(ell_theory[0])
    theory_ell_max = int(ell_theory[-1])
    data_ell_min = int(ell[0])
    data_ell_max = int(ell[-1])

    # Make sure we go up to high enough lmax
    if not (theory_ell_max >= data_ell_max):
        sys.stderr.write("You did not calculate the CMB to a high enough lmax needed lmax=%d but got %d\n" % (
            data_ell_max, theory_ell_max))
        return 1

    # Initialize to zero
    this_mod[:theory_ell_min, :] = 0.0

    # For correct BB window function we need theory BB (obviously) and theory EE
    #(lines 131-134 of bicep_util.py)
    # The "2" and "3" are for EE and BB respectively, to match the BICEP required ordering
    this_mod[theory_ell_min:data_ell_max + 1,
             2] = block[cmb_cl, "EE"][:data_ell_max - 1]
    this_mod[theory_ell_min:data_ell_max + 1,
             3] = block[cmb_cl, "BB"][:data_ell_max - 1]

    # generate the expected band-powers from the theory
    # using the window functions.
    expvBB = bu.calc_expvals(ell, this_mod, bpwf_l, bpwf_Cs_l)
    C_l[:, 0, 0] = expvBB[:, 3]

    # Add noise
    C_l = C_l + N_l

    # Get the likelihood
    like = float(bu.evaluateLikelihood(C_l, C_l_hat, C_fl, M_inv))

    # And save it
    block[likelihoods, "BICEP_LIKE"] = like

    # Signal that all is well
    return 0


def cleanup(config):
    # nothing to do here!  We just include this
    # for completeness
    return 0
