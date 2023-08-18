# COSEBIs

COSEBIs (Complete Orthogonal Sets of E/B-Integral) are a set of alternative two-point statistics 
designed to separate E/B-modes completely on a finite angular range. They have a reasonably localised
response to Fourier modes, ell, and are also easy to measure from data. This module calculates both $E_n$ and $B_n$ 
log-COSEBIs from [Schneider et al 2010](https://arxiv.org/abs/1002.2136). 

For this CosmoSIS edition of the COSEBIs library, we include only the cl_to_cosebis mode to enable KiDS-1000 or DES+KiDS
joint analyses.  This software is taken from Marika Asgari's [public COSEBIs library](https://github.com/maricool/2pt_stats).
Please visit Marika's library for additional CosmoSIS-compatible software that can, for example, convert finely binned $\xi_{\pm}(\theta)$ measurements into COSEBIs and calculate analytical covariance matrices.

The 'maricool' version of cl_to_cosebis has the capability to marginalise over an uncertain c-term.  We
do not include this functionality here for simplicity.

For B-mode enthusiasts, whilst the expectation value of B-modes is zero if there is no B-mode power spectra, there are B-modes from TATT intrinsic alignments.  You can calculate B-mode COSEBIs by simply switching the input_section_name to take B-mode Cls as input and output_section_name to avoid mixing E/B-modes.
