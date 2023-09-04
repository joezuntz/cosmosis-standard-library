#############
#############
Simon Samuroff
24th Aug 2023
#############
#############

This is a cosmosis wrapper for the BACCO emulator.
See this link for details about BACCO and links to the papers: https://baccoemu.readthedocs.io/en/latest/

The module itself uses the public baccoemu python package. I had problems getting the basic out-of-box pip installation to work, so I did the following

1. Source your cosmosis

2. Clone the repo
git clone https://bitbucket.org/rangulo/baccoemu.git

3. Use pip to install it from the repo
cd baccoemu
pip install . [--user]

After that it should be possible to import and use the emulator within a cosmosis pipeline

The module has two settings:

1. mode=pk_nl_only uses BACCO to get the nonlinear boost factor and applies it to a linear matter power spectrum already saved to the block

2. mode=pk_nl_plus_baryons uses BACCO to calculate both the nonlinear boost factor and the baryon suppression factor. Then it applies both to the linear P(k)

Notes:
The emulator has a limited range of k (0.001<k<5 Mpc/h) and z (<1.5). Beyond that we need to interpolate
