purpose: "implement baryonic effect in nonlinear pk based on hydrodynamic simulation measurements"
interface: baryonic_interface.py

example:
[eagle]
file = /usr/local/cosmosis/modules/baryonic_interface.py
mode = ratio 
; using hunjing's ratio fits file

ratiotable = logPkRatio_eagle.fits

[owls]
file = /usr/local/cosmosis/modules/baryonic_interface.py
mode = fixed
; using hunjing's ratio fits file

dmtable = powtable_DMONLY_all.dat
powtable = powtable_REF_all.dat

assumptions:
    - Pk scales linearly under baryonic effect. 
    - Pk_baryonic/Pk_nonlinear = Pk_hydro/Pk_DMonly

explanation: |
  this module reads in Pk_nonlinear from previous modules and output Pk_baryonic given hydro simulation pk+
  dmonly simulation data, or the ratio between them.
params:
    mode: "ratio" or "fixed"
    powtable: the hydrodynamic simulation pk in "fixed" mode.
    dmtable: the dmonly simulation pk in "fixed" mode.
    ratiotable: pk ratio fits file in "ratio" mode.
inputs:
    matter_power_nl:
        z: "1D real array, redshifts of samples"
        k_h: "1D real array, inpu k wavenumbers of samples in Mpc/h."
        p_k: "2D real array, matter power spectrum at samples in (Mpc/h)^-3."

outputs:
    matter_power_nl:
        z: "1D real array, redshifts of samples"
        k_h: "1D real array, inpu k wavenumbers of samples in Mpc/h, extended to kmax"
        p_k: "2D real array, matter power spectrum at samples in (Mpc/h)^-3, extended to kmax"
