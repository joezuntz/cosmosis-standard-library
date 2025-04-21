"""

Adapted from the module wigglez_bao and bao.f90 of cosmomc

"""
from cosmosis.datablock import option_section
import os
import numpy as np
from cosmosis.datablock import names as section_names

dist = section_names.distances
likes = section_names.likelihoods

ROOT_dir = os.path.split(os.path.abspath(__file__))[0]

dirname = os.path.split(__file__)[0]

# Depricated because now we have real DES Y6 BAO data
# MICE simulation (values from Juan Mena, computed on Jan 18 2023)
# If you want to use different values or use synthetic data vector created at a specific cosmology,\
# specify values at ini file. Otherwise we use the default quantities from MICE simulation
# REDSHIFT = 0.888 #(preliminary value)
# RS_FIDUCIAL = 153.4  #in Mpc
# DM_FIDUCIAL = 3110.6 #in Mpc

# Taken from DES Y6BAO likelihood 
# https://github.com/des-science/y6kp-bao-methods/blob/main/des-y6-bao/bao_y6_like.py
REDSHIFT = 0.851
# Fiducial Planck cosmology rs which was used to compute alpha
RS_FIDUCIAL =  147.8 
# Angular diameter distance at the fiducial cosmology, d_a(zeff=0.851)
DA_FIDUCIAL = 1627.70 
DM_FIDUCIAL = DA_FIDUCIAL * (1 + REDSHIFT)


def setup(options):
	section = option_section
	# file from here: https://github.com/des-science/y6kp-bao-methods/blob/main/des-y6-bao/chi2profile_dvdesy6_cosmoplanck18_covcosmolike.csv
	CHI2_file = os.path.join(ROOT_dir,'chi2profile_dvdesy6_cosmoplanck18_covcosmolike.txt')
	chi2 = np.loadtxt(CHI2_file)
	redshift = options.get_double(section, "redshift", default=REDSHIFT)
	rs_fiducial = options.get_double(section, "rs_fid", default=RS_FIDUCIAL)
	dm_fiducial = options.get_double(section, "dm_fid", default=DM_FIDUCIAL)
	feedback = options.get_int(option_section, "feedback", default=0)
	#Return data for later
	return (chi2,redshift,rs_fiducial,dm_fiducial,CHI2_file, feedback)

def execute(block, config):

	chi2,redshift,rs_fiducial, d_m_fiducial, CHI2_file, feedback = config
	
	alpha = chi2[:,0]
	if feedback:
		print("data_file = ", CHI2_file)
		print("alpha min = ", alpha[0])
		print("alpha max = ", alpha[-1])
	chi2_alpha = chi2[:,1] #careful with norm. Don't care about absolute value of likelihood . -0.5*np.log(2*np.pi*sigma**2))

	#load distance relations and R_S
	dist_z = block[dist, "z"]
	d_m = block[dist, 'd_m'] #in Mpc
	h = block[dist, "h"] #in 1/Mpc
	rs = block[dist, "RS_ZDRAG"] #in Mpc
	if dist_z[1] < dist_z[0]:
			dist_z = dist_z[::-1]
			d_m = d_m[::-1]
			h = h[::-1]


	#Compute the derived D_V distance
	#d_v = (dist_z*d_m*d_m/h)**(1./3.)

	#Interpolate the theory at the 
	#observed redshift
	d_m_predicted = np.interp(redshift, dist_z, d_m)
	alpha_predicted = (d_m_predicted/d_m_fiducial)*(rs_fiducial/rs)
	
	#If required, print out some info
	if feedback:
		print("redshift = ", redshift)
		print("alpha predicted = ",  alpha_predicted)
		print("rs = ", rs)
		print("rs_fiducial = ", rs_fiducial)
		print("dm = ", d_m_predicted)
		print("dm_fiducial = ", d_m_fiducial)

	#interpole the chi2 at the alpha predicted
	chi2_alpha_predicted = np.interp(alpha_predicted, alpha, chi2_alpha, left=chi2_alpha[0], right=chi2_alpha[-1])

	# if alpha_predicted < alpha[0]:
	# 	chi2_alpha_predicted = chi2_alpha[0]
	# elif alpha_predicted > alpha[-1]: 
	# 	chi2_alpha_predicted = chi2_alpha[-1]
	# else:
	# 	#simple interpolation of chi2 at the predicted alpha
	# 	ii = int(np.floor((alpha_predicted - alpha[0])/0.00404)) #0.00404 is the step of alpha in chi2 file
	# 	chi2_alpha_predicted = (chi2_alpha[ii+1] - chi2_alpha[ii])*(alpha_predicted-alpha[ii])/(alpha[ii+1] - alpha[ii]) + chi2_alpha[ii]

	#Get the log likelihood from the chi2
	like = -chi2_alpha_predicted/2.
	block[likes, 'DESY6BAO_LIKE'] = like

	#Signal success
	return 0

def cleanup(config):
	pass