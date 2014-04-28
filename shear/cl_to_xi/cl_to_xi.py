import numpy as np
import scipy.special
import scipy.integrate
import scipy.interpolate

def theta_range_for_ell(ell):
	#ranges 
	lmin = ell.min()
	lmax = ell.max()

	#Convert to output limits
	theta_min = np.pi/lmax
	theta_max = np.pi/lmin

	return theta_min, theta_max


def cl_to_xi_plus_and_minus(cl_interp,  cl_interp_type, ell_min, ell_max, thetas):
	# Make inputs into arrays and find limits we will work with
	#ell = np.array(ell)
	#c_ell = np.array(c_ell)
	lmin = ell_min#ell.min()
	lmax = ell_max#ell.max()
	#theta_min , theta_max = theta_range_for_ell(ell)

	# a few constant numbers.  We just increased these until things looked nice
	nzero = 20000
	# Set minimum number of points between zeros of integrand (will scale this number later)
	min_points_per_zero = 100
	#nell = ell.size

	# Zeros of the Bessel functions, with added 0 at the start as numpy does not 
	# We never actually use the full nzero range, we just keep going until 
	# we are past the regime where our c_ell are defined
	#include that one
	j0_zeros = scipy.special.jn_zeros(0, nzero)
	j4_zeros = scipy.special.jn_zeros(4, nzero)

	#Set up output space
	Theta = thetas #np.logspace(np.log10(theta_min), np.log10(theta_max), n_theta)
	Xi_plus = np.zeros_like(Theta)
	Xi_minus = np.zeros_like(Theta)

	#Do the integral separately for each theta value
	for theta_index, theta in enumerate(Theta):
		#zeros in ell
		ell_j0_zeros = j0_zeros/theta
		if lmin<ell_j0_zeros[0]:
			ell_j0_zeros=np.concatenate((np.array([lmin]),ell_j0_zeros))
		ell_j4_zeros = j4_zeros/theta
		if lmin<ell_j4_zeros[0]:
			ell_j4_zeros=np.concatenate((np.array([lmin]),ell_j4_zeros))
		xi_plus_theta = 0.0
		xi_minus_theta = 0.0
		#Get a fixed number of samples between each pair of zeros in the 
		#bessel function
		#Integrate xi_plus and xi_minus seperately as tey have different zeros
		#So first xi_plus:
		#Save integrand samples for debugging
		integrand_samples_debug=[]
		ells_debug=[]
		for i in xrange(nzero-1):
			ell_sample_min = ell_j0_zeros[i]
			ell_sample_max = ell_j0_zeros[i+1]
			# print ell_sample_min, ell_sample_max, lmin, lmax
			#We want to extrapolate at low l since Bessel functions
			#make integrand large here...don't think we need to extrapolate at high l
			#if there will never be any more then we are done
			#if ell_sample_max<lmin:
			#   continue
			if ell_sample_min>lmax:
				break
			if ell_sample_max>lmax:
				#This will be the last zero because next time
				#ell_sample_min will be the current ell_sample_max
				#so we should cut down the ell_sample_max to lmax.
				#We do not really need all the samples that we are 
				#using here but it will only hurt the runtime not
				#the accuracy.  This is important for the low theta
				#high ell regime
				ell_sample_max = lmax

			# sample points
			#print min_points_per_zero*(cl_interp(np.log10(0.5*(ell_sample_max+ell_sample_min))))/(-10)
			#points_per_zero=max(min_points_per_zero,min_points_per_zero*(cl_interp(np.log10(0.5*(ell_sample_max+ell_sample_min))))/(-10))
			ell_samples = ell_samples = np.linspace(ell_sample_min, ell_sample_max, min_points_per_zero)
			#print len(ell_samples)
			ells_debug+=list(ell_samples)
			#convert to x=ell*theta
			x_samples = ell_samples * theta
			#interpolate to get c_ell at these samples
			if cl_interp_type=='loglog':
				cl_samples=10**cl_interp(np.log10(ell_samples))
			elif cl_interp_type=='minus_loglog':
				cl_samples=-10**cl_interp(np.log10(ell_samples))
				#cl_samples[cl_samples>0.]=0.
			elif cl_interp_type=='log_ell':
				cl_samples=cl_interp(np.log10(ell_samples))
			else:
				print 'invalid cl_interp_type'
				return 1
			#lib function to get J0
			j0_samples = scipy.special.j0(x_samples)
			#Overall integrand
			integrand_samples = ell_samples * cl_samples * j0_samples
			integrand_samples_debug+=list(integrand_samples)
			#sample spacing
			dl = ell_samples[1] - ell_samples[0]
			# import pylab
			# pylab.plot(x_samples, j0_samples, ',')
			# pylab.plot(x_samples, cl_samples, ',')
			#Trapezium integration in this block
			delta = scipy.integrate.trapz(integrand_samples, dx=dl)
			#Add to cumulative total
			xi_plus_theta += delta
			#print delta
		#print ells_debug[:10]
		#print integrand_samples_debug[:10]
		#import pylab
		#pylab.plot(ells_debug,integrand_samples_debug,label='theta='+str(theta*(180./np.pi)*60))
		#pylab.legend()
		#pylab.xlabel('l')
		#pylab.ylabel('integrand')
		#pylab.show()
		#exit()
		#if theta_index==1:
		#       pylab.show()
		#   exit()
		integrand_samples_debug=[]
		ells_debug=[]
		for i in xrange(nzero-1):
			ell_sample_min = ell_j4_zeros[i]
			ell_sample_max = ell_j4_zeros[i+1]
			# print ell_sample_min, ell_sample_max, lmin, lmax
			#if there are no C_ell at this range then skip it.
			#if there will never be any more then we are done
			#if ell_sample_max<lmin:
			#   continue
			if ell_sample_min>lmax:
				break
			if ell_sample_max>lmax:
				#This will be the last zero because next time
				#ell_sample_min will be the current ell_sample_max
				#so we should cut down the ell_sample_max to lmax.
				#We do not really need all the samples that we are 
				#using here but it will only hurt the runtime not
				#the accuracy.  This is important for the low theta
				#high ell regime
				ell_sample_max = lmax

			# sample points
			ell_samples = np.linspace(ell_sample_min, ell_sample_max, min_points_per_zero)
			#convert to x=ell*theta
			x_samples = ell_samples * theta
			#interpolate to get c_ell at these samples
			if cl_interp_type=='loglog':
				cl_samples=10**cl_interp(np.log10(ell_samples))
			elif cl_interp_type=='minus_loglog':
				cl_samples=-10**cl_interp(np.log10(ell_samples))
			elif cl_interp_type=='log_ell':
				cl_samples=cl_interp(np.log10(ell_samples))
			else:
				print 'invalid cl_interp_type'
				return 1
			#lib function to get J0
			j4_samples = scipy.special.jn(4,x_samples)
			#Overall integrand
			integrand_samples = ell_samples * cl_samples * j4_samples
			#sample spacing
			dl = ell_samples[1] - ell_samples[0]
			# import pylab
			# pylab.plot(x_samples, j0_samples, ',')
			# pylab.plot(x_samples, cl_samples, ',')
			#Trapezium integration in this block
			delta = scipy.integrate.trapz(integrand_samples, dx=dl)
			#Add to cumulative total
			xi_minus_theta += delta

		# pylab.show()
		# x=raw_input()
		# if x.strip():
		#   sys.exit(0)
		#Summed value, with factor 2pi
		Xi_plus[theta_index] = xi_plus_theta / (2*np.pi)
		Xi_minus[theta_index] = xi_minus_theta / (2*np.pi)
	return np.array(Xi_plus), np.array(Xi_minus)

def cl_to_xi_plus_and_minus_extended_apodized(ell, c_ell, thetas):
	ell_min,ell_max=0.01,10000

	cl_interp,interp_type=get_cl_interp_function(ell,c_ell)
	
	#ell_sample=np.logspace(np.log10(ell_min),np.log10(ell_max),100)
	#pylab.loglog(ell_sample,10**cl_interp(np.log10(ell_sample)))
	#pylab.show()
	#exit()
	
	xi_plus, xi_minus = cl_to_xi_plus_and_minus(cl_interp, interp_type, ell_min, ell_max, thetas)
	return xi_plus, xi_minus

def get_cl_interp_function(ell,c_ell):
	"Get a interpolation function to use for cl"
	"if possible interpolate in loglog space"
	#print c_ell.min(),c_ell.max()
	if (c_ell>0).all():
		#print 'all c_ell > 0'
		cl_interp=extend_and_apodize(ell,c_ell)
		interp_type='loglog'
	elif (c_ell<=0).all():
		#print 'all c_ell < 0'
		interp_type='minus_loglog'
		cl_interp = extend_and_apodize(ell,-c_ell)
	else:
		interp_type='log_ell'
		cl_interp = scipy.interpolate.interp1d(np.log10(ell),c_ell,bounds_error=False,fill_value=0.)
	return cl_interp,interp_type    

def extend_and_apodize(ell, cl_bin, ell_min=0.01 , ell_max=10000, n_ell=100):
	"Extend a c_ell beyond its valid regime by fitting a 4-poly and apodize the extended bit"
	p_u=np.polyfit(np.log10(ell[-20:]), np.log10(cl_bin[-20:]), 2)
	p_l=np.polyfit(np.log10(ell[:10]), np.log10(cl_bin[:10]), 1)
	ell_sample = np.logspace(np.log10(ell_min), np.log10(ell_max), n_ell)
	cl_poly=np.zeros_like(ell_sample)
	#pylab.loglog(ell_sample,10**cl_poly)
	cl_poly[(ell_sample>=ell.min()) & (ell_sample<=ell.max())]=scipy.interpolate.interp1d(np.log10(ell),np.log10(cl_bin))(np.log10(ell_sample[(ell_sample>ell.min()) & (ell_sample<ell.max())]))
	cl_poly[ell_sample<ell.min()] = np.log10(cl_bin[0])#np.polyval(p_l, np.log10(ell_sample[ell_sample<ell.min()]))
	cl_poly[ell_sample>ell.max()] = np.polyval(p_u, np.log10(ell_sample[ell_sample>ell.max()]))
	cl_interp=scipy.interpolate.interp1d(np.log10(ell_sample),cl_poly,bounds_error=False,fill_value=0.)
	return cl_interp


def radians_to_arcmin(r):
	return np.degrees(r)*60

def arcmin_to_radians(a):
	return np.radians(a/60)

def test():
	import pyfits
	data = pyfits.getdata("./output.fits", "SHEAR_CL")
	ell = data.ELL
	c_ell = data.BIN6_1
	import pylab
	print 'start'
	# theta, xi = cl_to_xi_plus(ell, c_ell, len(ell))
	# pylab.loglog(theta, xi)
	# pylab.loglog(theta, -xi)
	theta_min, theta_max = theta_range_for_ell(ell)
	# print "Done 1"

	ell_ext, c_ell_ext = extend_and_apodize(ell, c_ell)
	theta_ext, xi_plus_ext, xi_minus_ext = cl_to_xi_plus_and_minus(ell_ext, c_ell_ext, len(ell_ext))
	print "Done"
	valid = (theta_ext>theta_min) & (theta_ext<theta_max)
	pylab.loglog(radians_to_arcmin(theta_ext[valid]), xi_plus_ext[valid],'-')
	pylab.loglog(radians_to_arcmin(theta_ext[valid]), -xi_plus_ext[valid],'-')
	pylab.loglog(radians_to_arcmin(theta_ext[valid]), xi_minus_ext[valid],'-')
	pylab.loglog(radians_to_arcmin(theta_ext[valid]), -xi_minus_ext[valid],'-')

	pylab.show()
# 

# import cProfile
# cProfile.run('test()', sort='tottime')
# test()

