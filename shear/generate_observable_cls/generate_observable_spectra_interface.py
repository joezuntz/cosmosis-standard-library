from __future__ import print_function
from cosmosis.datablock import names, option_section
import generate_observable_spectra as functions

import numpy as np
import matplotlib.pyplot as plt


def setup(options):
    # DO options for theory terms
    shear = options[option_section, "shear"]
    intrinsic_alignments = options[option_section, "intrinsic_alignments"]
    cmb_kappa = options[option_section, "cmb_kappa"]
    kappa_shear = options[option_section, "kappa_shear"]
    kappa_pos = options[option_section, "kappa_position"]
    if intrinsic_alignments:
        GI = options.get_bool(option_section, "GI", default=True)
        II = options.get_bool(option_section, "II", default=True)
    else:
        GI = False
        II = False
    position = options[option_section, "position"]
    ggl = options[option_section, "ggl"]
    magnification = options[option_section, "magnification"]

    # DO options for binning and shot noise
    noise = options.get_bool(option_section, "noise", default=False)
    bias = (options[option_section, "bias"],
            options[option_section, "m_per_bin"])
    bins = options[option_section, "angular_frequency_bins"]
    window = options.get_string(option_section, "window", default="")

    nlbin_shear = options.get_int(option_section, "nlbin_shear", default=0)
    lmax_sh = options.get_double(option_section, "lmax_shear", default=0.)
    lmin_sh = options.get_double(option_section, "lmin_shear", default=0.)

    nlbin_ggl = options.get_int(option_section, "nlbin_ggl", default=0)
    lmax_ggl = options.get_double(option_section, "lmax_ggl", default=0.)
    lmin_ggl = options.get_double(option_section, "lmin_ggl", default=0.)

    nlbin_pos = options.get_int(option_section, "nlbin_pos", default=0)
    lmax_pos = options.get_double(option_section, "lmax_pos", default=0.)
    lmin_pos = options.get_double(option_section, "lmin_pos", default=0.)

    nlbin_kk = options.get_int(option_section, "nlbin_cmb_kappa", default=0)
    lmax_kk = options.get_double(option_section, "lmax_cmb_kappa", default=0.)
    lmin_kk = options.get_double(option_section, "lmin_cmb_kappa", default=0.)

    nlbin_ke = options.get_int(option_section, "nlbin_kappa_shear", default=0)
    lmax_ke = options.get_double(
        option_section, "lmax_kappa_shear", default=0.)
    lmin_ke = options.get_double(
        option_section, "lmin_kappa_shear", default=0.)

    nlbin_kn = options.get_int(
        option_section, "nlbin_kappa_position", default=0)
    lmax_kn = options.get_double(
        option_section, "lmax_kappa_position", default=0.)
    lmin_kn = options.get_double(
        option_section, "lmin_kappa_position", default=0.)

    disp = {True: 'yes', False: 'no'}
    print('Shot noise: %s' % disp[noise])
    print('Shape measurement bias (in each bin): %s (%s)' % (disp[bias[0]], disp[bias[1]]))
    print('Angular frequency binning: %s' % disp[bins])
    if bins:
        print('Using window function: %s' % window)

    shear_survey = options.get_string(
        option_section, "shear_sample", default="")
    pos_survey = options.get_string(option_section, "lss_sample", default="")
    cmb_survey = options.get_string(option_section, "cmb_sample", default="")
    if shear_survey:
        ngal_shear = options.get_double(shear_survey, "ngal", default=0.)
    if pos_survey:
        ngal_pos = options.get_double(pos_survey, "ngal", default=0.)
    sigma_gamma = options.get_double(
        shear_survey, "shape_dispersion", default=0.25)

    output_datavector = options.get_string(
        option_section, "output", default="")

    opt = {'shear': shear, 'nlbin_shear': nlbin_shear, 'nlbin_ggl': nlbin_ggl, 'nlbin_pos': nlbin_pos,
           'lmax_shear': lmax_sh, 'lmin_shear': lmin_sh,
           'lmax_pos': lmax_pos, 'lmin_pos': lmin_pos,
           'lmax_ggl': lmax_ggl, 'lmin_ggl': lmin_ggl,
           'lmax_kk': lmax_kk, 'lmin_kk': lmin_kk,
           'lmax_ke': lmax_ke, 'lmin_ke': lmin_ke,
           'lmax_kn': lmax_kn, 'lmin_kn': lmin_kn,
           'intrinsic_alignments': intrinsic_alignments, 'GI': GI, 'II': II,
           'position': position,
           'ggl': ggl,
           'magnification': magnification,
           'cmb_kappa': cmb_kappa,
           'kappa_shear': kappa_shear,
           'kappa_pos': kappa_pos,
           'noise': noise, 'bias': bias,
           'binning': bins, "window": window,
           'shear_cat': shear_survey, 'pos_cat': pos_survey, 'cmb_survey': cmb_survey,
           'output_datavector': output_datavector,
           'ngal': (ngal_shear, ngal_pos),
           'shape_dispersion': sigma_gamma}
    return opt


def execute(block, config):

    out_path = config['output_datavector']

    Cl = functions.Cl_class(block, config)
    Cl.load_and_generate_observable_cls(block, names)
    Cl.save_cls(block, out_path)

    return 0


def cleanup(config):
    pass
